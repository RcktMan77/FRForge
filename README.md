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

| Status | Milestone / phase |
|--------|-------------------|
| **Completed** | **M0–M8** — 1D/2D FR laboratory, invent/score, HO VTK |
| **In progress** | **Phase 2** — Research infrastructure & robustness |
| **Completed** | **Phase 2** (experiment log · schemes · robustness matrix) |
| **In progress** | **Phase 3** — Scientific & geometric extensions |
| **Completed** | **P3.1** Richer 2D capturing |
| **Completed** | **P3.2** Curved elements |
| **Completed** | **P3.3a** 2D Riemann + isentropic vortex |
| **Completed** | **P3.3b** Optional reduced DMR/FFS |
| **Completed** | **P4** Performance (alloc reduction, numeric fidelity) |
| **In progress** | **P5.1** Experiment-log analytics |
| Next | **P5.2** Publication / reproducibility snapshots |

### Milestone roadmap (M0–M8 complete)

| # | Name | Success criterion (summary) |
|---|------|------------------------------|
| 0 | Skeleton & harness | `frforge test --report` → valid JSON; green CI |
| 1 | 1D FR linear advection | Order ≈ formal (within 0.3) for \(p=2,3,4\); periodic mass ~ machine precision |
| 2 | 1D inviscid Burgers | Runs, conserves, shows HO oscillations |
| 3 | 1D Euler + smooth order | Formal order with capturing off; non-periodic BC path |
| 4 | Pluggable capturing interface | Abstract hooks + Persson AV baseline |
| 5 | Quantitative suite | Sod, Shu–Osher (self-converged hashed ref); scored JSON |
| 6 | Invention loop | `invent` / `score` vs baseline; `promising` / `accepted_candidate` status |
| 7 | High-order VTK writer | ParaView-ready HO VTU (may start after M3 for debug) |
| 8 | 2D + visualization | 2D order tests + HO VTK |

### Phase 2+ (post-M8)

| ID | Name | Summary |
|----|------|---------|
| P2.1 | Experiment log | `research/experiment_log.md` + invent auto-append |
| P2.2 | Configurable schemes | GLL, HLLC, SSP-RK2; invent defaults frozen |
| P2.3 | Robustness matrix | `frforge robustness --matrix ci\|full`; promotion rule |
| P3.x | Science/geometry | 2D capturing, curved elements, extra benchmarks |
| P4 | Performance | Profile / threads after research workflow is stable |

**Agent invent workflow:** (1) read [`research/experiment_log.md`](research/experiment_log.md) and/or `frforge log summary|frontier|lessons`; (2) implement under `src/methods/`; (3) `frforge invent --method …`; (4) fill hypothesis/lessons if promising+; (5) when short-listed, freeze a reproducibility snapshot (P5.2).

---

## Requirements

- Julia **≥ 1.10** (developed/tested on 1.11.x)
- Git

### CI & platforms

| Environment | Role |
|-------------|------|
| **Ubuntu (GitHub Actions)** | **Primary / required / blocking** CI on PRs to `develop` and `main` (Julia 1.10 and 1.11) |
| **macOS (local)** | **First-class** for development and verification — the dependency set (stdlib, JSON, ArgParse) is highly portable |
| **macOS (GHA)** | Optional, non-blocking job later once core milestones are stable — not required for early milestones |

**CI two-tier policy:** Every addition declares **required CI** vs **full/nightly/manual**. Required CI = default scheme only, modest grids, reduced Shu–Osher; **no** full robustness matrix or heavy 2D benchmarks. Budget: **~10–15 min** on Ubuntu. No large VTU/invent trees as PR artifacts.

---

## Setup

```bash
git clone https://github.com/RcktMan77/FRForge.git
cd FRForge
julia --project=. -e 'using Pkg; Pkg.instantiate()'
chmod +x bin/frforge
```

Works on macOS and Linux with the same commands.

---

## CLI

```bash
# Smoke: valid JSON skeleton only
./bin/frforge test --suite smoke --report results/smoke/report.json

# Milestone 1: linear advection order + conservation (p=2,3,4)
./bin/frforge test --suite advection --report results/m1/report.json

# Milestone 2: Burgers conservation + HO oscillation demo
./bin/frforge test --suite burgers --report results/m2/report.json

# Milestone 3: Euler density-wave order + BC freestream tests
./bin/frforge test --suite euler --report results/m3/report.json

# Milestone 4: Persson AV baseline vs NullCapturing on Burgers
./bin/frforge test --suite capturing --report results/m4/report.json

# Milestone 5: quantitative suite (order + Sod + Shu–Osher) with scored JSON
./bin/frforge test --suite quant --method persson_av --report results/m5/report.json

# Milestone 6: invent / score a method vs Persson baseline
./bin/frforge invent --method scaled_persson --baseline persson_av --report-dir results/invent

# Milestone 7–8: high-order VTU (1D and 2D) for ParaView (≥ 5.5)
./bin/frforge run --case euler_density_wave --p 3 --ne 16 --output results/euler.vtu
./bin/frforge run --case euler2d_wave --p 2 --ne 8 --t-final 0.5 --output results/euler2d.vtu
./bin/frforge test --suite 2d --report results/m8/report.json

# Help
./bin/frforge --help
```

| Command | Purpose | Available |
|---------|---------|-----------|
| `frforge test [--report PATH] [--suite …] [--method …]` | Verification → scored JSON | M0+ |
| `frforge run --case … [--method …]` | Single run | M1+ |
| `frforge invent --method M --baseline persson_av` | Quant suite + `candidate_status` | M6 |
| `frforge score --method-report a.json --baseline-report b.json` | Classify two reports | M6 |
| `frforge run … --output sol.vtu` | High-order discontinuous Lagrange VTU | M7 |

### ParaView notes

- Prefer **ParaView ≥ 5.5** for `VTK_LAGRANGE_LINE` (type 68) support.
- 1D solutions are written as discontinuous high-order line cells (jumps visible at faces).
- Euler fields: `rho`, `u`, `p` (primitives) and `rho`, `rho_u`, `E` (conserved).

---

## How to add a new method

1. **Create** `src/methods/MyMethod.jl` implementing `AbstractCapturingMethod` (override only the hooks you need: `preprocess_state!`, `extrapolate_interface!`, `numerical_flux_method`, `sense!`, `apply_dissipation!`, `post_step!`).
2. **Include + register** in `src/methods/Registry.jl`:
   ```julia
   include("MyMethod.jl")
   register_method!("my_method", (; kwargs...) -> MyMethod(; kwargs...))
   ```
3. **Run** the invent loop:
   ```bash
   ./bin/frforge invent --method my_method --baseline persson_av --report-dir results/invent
   ```
4. **Read** `candidate_status` in the printed summary / JSON:
   - `rejected` — hard gates failed  
   - `pass_gates` — valid but not better than baseline by δ (default 0.02) or trade-off failed  
   - `promising` — better composite + order-vs-dissipation trade-off  
   - `accepted_candidate` — promising **and** high-order VTK produced (M7+)  
5. Structural novelty (new sensors/operators/hybrids) is the primary research goal; coefficient-only search is allowed but secondary.

See also `docs/design.md` § Agent invention pathway.

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
