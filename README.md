# FRForge

**High-Order Flux Reconstruction Laboratory for Novel Shock-Capturing Methods**

FRForge is a green-field, self-contained [Julia](https://julialang.org/) package for inventing and rigorously evaluating new high-order shock-capturing / discontinuity-treatment methods for Flux Reconstruction (FR) schemes applied to the compressible Euler equations.

The laboratory is built so that agents (and humans) can propose **structurally novel** methods — new sensors, dissipation operators, hybrid schemes, limiters — not merely tune coefficients of a fixed heuristic. Every candidate is judged by a quantitative suite that measures:

| Metric | Intent |
|--------|--------|
| Formal order of accuracy | Smooth-region fidelity |
| Excess numerical dissipation | Do not over-smooth |
| Shock resolution & overshoots | Capture quality |
| Positivity & conservation | Robustness / correctness |

Full design blueprint: [`docs/design.md`](docs/design.md).

---

## Invention philosophy

1. **Structural novelty first** — new code against abstract residual hooks, not only parameter search.
2. **Machine-readable verdicts** — every verification run emits JSON agents can score autonomously.
3. **Baseline honesty** — a classical Persson-style modal sensor + artificial viscosity is a reference, not the goal.
4. **Milestone discipline** — each milestone ends with runnable code and passing verification before the next begins.
5. **Clarity over performance** (early) — pure Julia, first-principles FR, readable interfaces.

---

## Current milestone

| Status | Milestone |
|--------|-----------|
| **Completed** | **0 — Repository setup, project skeleton & verification harness** |
| **Next** | **1 — 1D FR linear advection** |

### Milestone roadmap

| # | Name | Success criterion (summary) |
|---|------|------------------------------|
| 0 | Skeleton & harness | `frforge test --report` → valid JSON; green CI |
| 1 | 1D FR linear advection | Order ≈ formal (within 0.3); tight conservation |
| 2 | 1D inviscid Burgers | Runs, conserves, shows HO oscillations |
| 3 | 1D Euler + smooth order | Formal order with capturing off |
| 4 | Pluggable capturing interface | Abstract hooks + Persson AV baseline |
| 5 | Quantitative suite | Sod, Shu–Osher, scored JSON summary |
| 6 | Invention loop | `invent` / `score` vs baseline |
| 7 | High-order VTK writer | ParaView-ready HO VTU |
| 8 | 2D + visualization | 2D order tests + HO VTK |

---

## Requirements

- Julia **≥ 1.10** (developed/tested on 1.11.x)
- Git

---

## Setup

```bash
git clone https://github.com/RcktMan77/FRForge.git
cd FRForge
julia --project=. -e 'using Pkg; Pkg.instantiate()'
chmod +x bin/frforge
```

---

## CLI

```bash
# Milestone 0: write a valid JSON report skeleton
./bin/frforge test --report results/smoke/report.json

# Help
./bin/frforge --help
```

| Command | Purpose | Available |
|---------|---------|-----------|
| `frforge test [--report PATH] [--suite NAME]` | Verification → JSON | M0+ |
| `frforge run ...` | Single case (+ optional VTK) | M1+ / M7 |
| `frforge invent ...` | Method vs baseline | M6 |
| `frforge score ...` | Score two reports | M6 |

Reports land under `results/` by convention (gitignored). Schema is versioned (`schema_version: 1`); see design doc for field tables and scoring weights.

---

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Or via the CLI smoke path above (also exercised in CI).

---

## Git workflow

- **`main`** — protected, always green/releasable. Never commit directly.
- **`develop`** — integration branch.
- **`milestone/*` / `feature/*`** — short-lived work branches; prefer PRs (even self-merged).

```text
milestone/N-...  →  develop  →  main
```

---

## License

[MIT](LICENSE) © 2026 RcktMan77
