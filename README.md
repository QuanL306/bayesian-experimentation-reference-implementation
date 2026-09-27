# Bayesian Experimentation: Reference Implementation

End-to-end Bayesian A/B testing, multi-armed bandit simulation, and experiment design comparison — implemented from scratch in R with publication-quality visualizations.

## Overview

This repository contains three self-contained projects that implement the core methodologies used by production experimentation platforms (Eppo, Statsig, GrowthBook, Netflix Ax). Every algorithm is implemented from first principles — no wrapper libraries.

| Project | What It Does | Key Methods |
|---------|-------------|-------------|
| [01-bayesian-ab-testing](01-bayesian-ab-testing/) | Full Bayesian A/B test workflow for e-commerce checkout | Beta-Binomial conjugate model, ROPE, sequential monitoring, prior sensitivity |
| [02-multi-armed-bandits](02-multi-armed-bandits/) | 8 algorithms × 5 scenarios systematic benchmark | Thompson Sampling, UCB variants, EXP3, non-stationary adaptations |
| [03-experiment-design](03-experiment-design/) | Power & efficiency comparison across experiment paradigms | Bayesian vs frequentist power, multiple comparisons, switchback design |

## 01 · Bayesian A/B Testing Pipeline

A complete decision workflow for an e-commerce checkout optimization (control 8% → treatment 9% conversion).

**Sections:**
1. **Pre-experiment power analysis** — simulation-based Bayesian power under ROPE framework
2. **Posterior inference** — Beta-Binomial conjugate updates with 95% HDI
3. **ROPE decision** — three-way decision (accept H₀ / reject H₀ / uncertain) with δ = ±0.5pp
4. **Sequential monitoring** — daily posterior updates with O'Brien-Fleming-style boundaries
5. **Subgroup analysis** — forest plot of treatment effects across 5 user segments
6. **Prior sensitivity** — 5 prior specifications from non-informative to informative
7. **Revenue posterior** — joint modeling of conversion + ARPU for business impact

> **7 figures** · ~500 lines R · `notes.md` with 7 sections of methodology documentation

**Run:** `Rscript 01-bayesian-ab-testing/bayesian_ab_fullflow.R`

---

## 02 · Multi-Armed Bandit Benchmark

Systematic comparison of 8 bandit algorithms across 5 difficulty scenarios.

**Algorithms (all from scratch):**
- ε-greedy, UCB1, UCB-Tuned, Thompson Sampling, EXP3
- Sliding-Window TS, Discounted UCB, UCB-V (variance-aware)

**Scenarios:**
| Scenario | Challenge | Best Algorithm |
|----------|-----------|---------------|
| Standard (K=5, clear gaps) | Baseline | Thompson Sampling |
| Hard (K=5, tiny gaps) | Distinguishing similar arms | UCB-Tuned |
| Many arms (K=20) | Exploration cost | Thompson Sampling |
| Non-stationary (drift) | Adapting to change | Sliding-Window TS |
| Gaussian rewards | Reward distribution | Thompson Sampling |

**Key finding:** Thompson Sampling achieves sublinear regret O(Σ log T / Δᵢ) and outperforms UCB1 by 15–30% in finite-sample cumulative regret across most scenarios. Non-stationary variants (SW-TS, Discounted UCB) recover quickly after regime shifts where static algorithms fail.

> **10 figures** · ~600 lines R · `notes.md` with algorithm derivations and regret analysis

**Run:** `Rscript 02-multi-armed-bandits/bandit_simulation.R`

---

## 03 · Experiment Design Comparison

When should you use a standard A/B test vs. multi-variant vs. switchback design?

**Three scenarios analyzed:**

**Scenario 1 — E-commerce A/B Test**
- Frequentist analytical power vs. Bayesian simulation power
- Sample size requirements across MDEs (0.5pp – 3pp)
- Bayesian power under ROPE is more conservative than frequentist power at same n

**Scenario 2 — Multi-Variant (4 arms)**
- Bonferroni, Holm-Bonferroni, and Bayesian hierarchical correction
- FWER inflation: 20 tests at α=0.05 → 64% chance of ≥1 false positive
- Bayesian hierarchical model naturally shrinks extreme estimates

**Scenario 3 — Switchback Experiment**
- Time-block randomization for marketplace/platform interventions
- Carryover effect detection (lag-1 through lag-5)
- Naive SE vs. HAC SE: naive underestimates by ~1.2×
- Optimal block length trades off carryover vs. statistical efficiency

> **5 figures** · ~400 lines R · `notes.md` with design decision framework

**Run:** `Rscript 03-experiment-design/experiment_design.R`

---

## Methodology Decisions

These projects make deliberate methodological choices that reflect production experimentation practice:

- **ROPE over p-values**: Decisions incorporate effect-size relevance (δ chosen by business context), not just statistical significance
- **Sequential monitoring**: Bayesian framework supports optional stopping without α-inflation, with O'Brien-Fleming boundaries to prevent premature decisions
- **Thompson Sampling over UCB**: Better finite-sample performance, natural prior incorporation, probabilistic outputs easier to communicate to stakeholders
- **HAC standard errors for switchback**: Accounts for temporal autocorrelation; naive SE underestimates uncertainty by ~20%
- **Bayesian hierarchical models for multiple comparisons**: Avoids the conservativeness of Bonferroni while naturally controlling false discoveries through shrinkage

## Tech Stack

- **Language:** R 4.x
- **Dependencies:** ggplot2, gridExtra (visualization only — all statistics implemented from scratch)
- **No proprietary data:** All data is simulated with known ground truth, enabling exact verification of method performance

## How to Run

```bash
# Clone
git clone https://github.com/QuanL306/bayesian-experimentation-reference-implementation.git
cd bayesian-experimentation-reference-implementation

# Run any project
Rscript 01-bayesian-ab-testing/bayesian_ab_fullflow.R
Rscript 02-multi-armed-bandits/bandit_simulation.R
Rscript 03-experiment-design/experiment_design.R
```

Each script is self-contained: generates its own simulated data, runs all analyses, and saves figures to its own directory.

## Related Repositories

| Repository | Focus |
|-----------|-------|
| [`causal-inference-python`](https://github.com/QuanL306/causal-inference-python) | Python causal inference: PS/IPW/DR, DML, Causal Forests, DID, RDD, IV, causal discovery |
| [`multi-method-causal-replication`](https://github.com/QuanL306/multi-method-causal-replication) | Signature capstone: end-to-end causal pipeline with 6-method triangulation |

## License

MIT
