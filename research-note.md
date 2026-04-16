# Klima Protocol — Fee-Based vs. Principal Model

A research note on revenue and risk quantification for a carbon-retirement intermediary operating against the Klima Protocol.

This is **Phase A** of a three-phase roadmap:

- Phase A (this note) — mathematical framework, closed-form moments under GBM, risk-metric menu.
- Phase B — TypeScript simulation reproducing the results below.
- Phase C — jump-diffusion / regime-switching extensions, Poisson demand, historical calibration against a kVCM proxy, hedging strategies.

## 1. Setup and notation

| Symbol | Meaning |
| --- | --- |
| `S_t` | kVCM/USD spot price at time `t` |
| `P` | Protocol price (constant) in kVCM per tonne |
| `π_t := P · S_t` | USD protocol price per tonne |
| `λ` | Deterministic retirement flow, tonnes / unit time |
| `T` | Horizon |
| `N := λ · T` | Total tonnes retired over `[0, T]` |
| `f` | Fee rate (e.g. 0.05) |
| `Q` | Fixed USD quote per tonne in the principal model |
| `α ∈ [0, 1]` | Fraction of inventory pre-purchased (principal, generalized) |

**Price dynamics.** `S_t` follows Geometric Brownian Motion under the physical measure:

```
dS_t = μ · S_t · dt + σ · S_t · dW_t,    S_0 given.
```

Because `P` is a constant multiplier, `π_t = P · S_t` is itself GBM with the same `(μ, σ)`.

**Demand.** Retirement flow is deterministic at rate `λ` tonnes / unit time. The sole source of randomness in the baseline is the kVCM price.

**Ignored in the baseline.** Risk-free discounting, gas, on-chain slippage, order-flow stochasticity. Each is re-introduced in a later phase (see §6).

### The central stochastic object

Both models' P&L are linear functionals of the **integral of the GBM**:

```
I_T := ∫₀ᵀ S_t dt.
```

Its first two moments are known in closed form (Dufresne, 2001):

```
E[I_T]  = S_0 · (e^{μT} − 1) / μ,                                                   μ ≠ 0    (→ S_0 · T as μ → 0)

E[I_T²] = (2 · S_0²) / (μ + σ²)
          · [ (e^{(2μ + σ²)T} − 1) / (2μ + σ²)  −  (e^{μT} − 1) / μ ],

Var[I_T] = E[I_T²] − E[I_T]².
```

The distribution of `I_T` is not log-normal, so tail metrics (VaR, CVaR) are obtained by Monte Carlo; the moments above serve as closed-form anchors and Monte Carlo sanity checks.

## 2. Fee-based model

The company quotes clients `(1 + f) · π_t` and remits `π_t` to the protocol, keeping `f · π_t` per tonne.

**Total revenue over `[0, T]`:**

```
R_fee = ∫₀ᵀ f · π_t · λ dt = f · P · λ · I_T.
```

**Moments:**

```
E[R_fee]   = f · P · λ · E[I_T].
Var[R_fee] = (f · P · λ)² · Var[I_T].
```

**Properties.**

- `R_fee ≥ 0` almost surely.
- Top-line volatility is fully driven by `σ` — the company has no balance-sheet exposure and holds no inventory.
- Scales linearly in `f`, so the risk-adjusted return per unit fee is invariant in `f`.

## 3. Principal model

The company sets a fixed USD quote `Q` at `t = 0`. Three variants of inventory sourcing, in increasing order of risk:

### 3a. Fully matched pre-purchase (α = 1)

At `t = 0` buy exactly `N · P` kVCM at spot `S_0`; cost `C = N · P · S_0`. Burn against deterministic demand over `[0, T]`.

**Terminal P&L:**

```
Π_matched = N · (Q − P · S_0),
```

which is **deterministic** — no terminal kVCM risk.

Risk still exists *interim*. The mark-to-market inventory value at time `t` is

```
V_t = (N · P) · (1 − t/T) · S_t,
```

which is a scaled GBM decayed by a deterministic burn schedule. Solvency, margin-call, or accounting-covenant concerns live in the distribution of `max_{t ≤ T} (V_0 − V_t)` (max drawdown). This is tracked by Monte Carlo in Phase B.

### 3b. Back-to-back acquisition (α = 0)

The company quotes `Q` at `t = 0` but buys kVCM at spot for each retirement. Per-tonne realized P&L is `Q − P · S_t`.

**Total P&L:**

```
Π_b2b = ∫₀ᵀ (Q − P · S_t) · λ dt = Q · N − P · λ · I_T.
```

**Moments:**

```
E[Π_b2b]   = Q · N − P · λ · E[I_T].
Var[Π_b2b] = (P · λ)² · Var[I_T].
```

Note that `Var[Π_b2b]` coincides with `Var[R_fee]` up to the rescaling factor `(P / f)² · f² = P²` — i.e. the two models share **the same random kernel** `I_T`. Phase B can therefore reuse a single Monte Carlo simulation of `I_T` and rescale.

**Payoff shape.** `Π_b2b` is linearly *decreasing* in `I_T`: upside is capped at `Q · N` (reached as `S_t → 0`), downside is unbounded if kVCM rallies. Equivalent to shorting a continuous strip of forwards on kVCM struck at `Q / P`.

### 3c. Partial pre-purchase (α ∈ [0, 1])

Buy `α · N · P` kVCM at `S_0`, source the rest back-to-back. Then

```
Π_α = α · N · (Q − P · S_0)  +  (1 − α) · (Q · N − P · λ · I_T)
     = Q · N  −  P · λ · [ α · S_0 · T  +  (1 − α) · I_T ].
```

Mean and variance interpolate linearly in `α`:

```
E[Π_α]   = (1 − α) · E[Π_b2b] + α · Π_matched,
Var[Π_α] = (1 − α)² · Var[Π_b2b].
```

`α` is the company's hedge ratio against spot. `α = 1` gives a deterministic P&L (fully hedged at inception); `α = 0` is the unhedged short-forward strip.

## 4. Risk quantification

For each model, the Phase B simulator should report:

| Metric | How |
| --- | --- |
| `E[Π]`, `Var[Π]`, `SD[Π]` | Closed form from §2–3 |
| `VaR_95`, `VaR_99` | Monte Carlo empirical quantile of `−Π` |
| `CVaR_95`, `CVaR_99` | Monte Carlo tail mean of `−Π` |
| `P[Π < 0]` | Monte Carlo |
| `Sharpe-like = E[Π] / SD[Π]` | Closed form |
| Max NAV drawdown (principal 3a only) | Monte Carlo on `V_t` path |

### Itô dynamics and delta

Applying the Itô product rule to the cumulative P&L process, the instantaneous sensitivity of *remaining* P&L to the spot `S_t` is:

```
Fee-based:     ∂ E[R_fee − R(t) | F_t] / ∂ S_t  =  f · P · λ · (e^{μ(T−t)} − 1) / μ
                                                 ≈ f · P · λ · (T − t)   for μT small.

Principal 3b:  ∂ E[Π_b2b − Π(t) | F_t] / ∂ S_t  =  − P · λ · (e^{μ(T−t)} − 1) / μ
                                                 ≈ − P · λ · (T − t).
```

Signs are opposite: the fee book is **long** kVCM beta; the back-to-back principal book is **short** kVCM beta. The matched principal book has zero delta (fully pre-hedged by physical inventory).

This is the handle for Phase C hedging: the natural static hedge for the principal back-to-back book is to hold `(P · λ) · (T − t)` tokens of spot kVCM at each time `t` — which is exactly the matched-pre-purchase strategy (§3a) amortized to the remaining horizon.

## 5. Direct comparison

| | Fee-based | Principal 3a (matched) | Principal 3b (back-to-back) |
| --- | --- | --- | --- |
| `E[P&L]` | `f · P · λ · E[I_T]` | `N · (Q − P · S_0)` | `Q · N − P · λ · E[I_T]` |
| `Var[P&L]` | `(f P λ)² · Var[I_T]` | `0` (terminal) | `(P λ)² · Var[I_T]` |
| kVCM exposure | long | none (terminal), long (interim NAV) | short |
| Downside | bounded below by 0 | deterministic | unbounded |
| Capital requirement | none | `N · P · S_0` | none |

### Break-even quote

The principal back-to-back model matches the fee-based model's *expected* revenue when

```
Q* = (1 + f) · P · E[I_T] / T = (1 + f) · P · S_0 · (e^{μT} − 1) / (μ T).
```

As `μ → 0`, `Q* → (1 + f) · P · S_0` — the fee-based time-zero quote. For `μ > 0` the principal model must quote *above* that to compensate for expected kVCM appreciation; for `μ < 0` it quotes below.

**Asymmetry observation.** Even at `Q = Q*`, the two books have radically different risk profiles: the fee book has bounded positive-only revenue, while the back-to-back principal book has the same variance but a *left-skewed* loss tail. Matching means does not match distributions.

## 6. Limitations and next steps

| Baseline simplification | Removed in |
| --- | --- |
| Deterministic demand | Phase C — compound-Poisson order flow |
| GBM price dynamics | Phase C — Merton jump-diffusion and/or two-state regime switching |
| No historical calibration | Phase C — kVCM proxy (KLIMA, BCT, NCT) |
| No discounting, gas, or on-chain slippage | Phase B — parameterized |
| Static (or absent) hedging | Phase C — dynamic delta hedge with inventory; perp/futures hedge if available |
| No credit / counterparty layer | Not scoped |

## References

- Dufresne, D. (2001). *The integral of geometric Brownian motion.* Advances in Applied Probability, 33(1), 223–241. — closed-form moments of `I_T`.
- Glasserman, P. (2003). *Monte Carlo Methods in Financial Engineering*, §3.4. — simulation of path integrals of GBM.
