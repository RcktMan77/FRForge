# FRForge

**A laboratory for inventing and quantitatively evaluating high-order shock-capturing methods in Flux Reconstruction.**

FRForge is a self-contained [Julia](https://julialang.org/) package for the compressible Euler equations that treats *discontinuity treatment* as the research object: not only whether a scheme runs, but whether a proposed sensor, artificial viscosity, hybrid reconstruction, or limiter improves on a classical baseline without destroying formal accuracy in smooth flow. It is intended for researchers and numerical analysts who want a transparent, first-principles FR core, a pluggable capturing interface, and machine-readable scores that make method comparison reproducible.

---

## The problem

High-order discontinuous schemes—Flux Reconstruction (FR), discontinuous Galerkin (DG), and spectral difference methods—deliver excellent accuracy and dispersion properties when the solution is smooth. At shocks and contacts, however, polynomial representations produce Gibbs-type oscillations that can destroy positivity of density and pressure and crash the computation.

Classical remedies include slope/WENO limiting, filtering, and **element-local artificial viscosity** guided by a modal sensor (Persson & Peraire and related work). These approaches are successful in practice, but leave well-known trade-offs:

- Viscosity strong enough for robustness often **over-dissipates** entropy waves and fine post-shock structure.
- Sensors and viscosity scales that work on one mesh or polynomial degree may fail on another.
- Hybrid or nonlinear reconstructions can restore sharpness, yet may compromise **formal order** in smooth regions if not carefully designed.
- Literature comparisons are hard to reproduce: different codes, mesh families, time integrators, and qualitative “looks good” criteria.

What is still scarce is a **systematic laboratory**: a fixed high-order FR discretization, a clear capturing hook surface for *structurally new* methods (not only coefficient tweaks), and a quantitative suite that scores order preservation, excess dissipation, shock quality, and robustness against a documented baseline—so that candidate methods can be accepted, rejected, or short-listed with evidence. The concrete goal is to **produce methods that improve the measured order–dissipation–shock trade-off relative to a classical Persson AV baseline under a fixed FR discretization**.

---

## Approach

FRForge is organized around four ideas:

1. **Green-field high-order FR core**  
   Solution points, correction functions recovering a DG-equivalent operator, interface numerical fluxes, and strong-stability-preserving Runge–Kutta time integration are implemented from first principles in pure Julia (readable operators, no opaque third-party FR dependency).

2. **Pluggable capturing interface**  
   Capturing methods implement abstract residual hooks (`preprocess_state!`, interface extrapolation, optional flux override, `sense!`, `apply_dissipation!`, `post_step!`). New sensors, dissipation operators, hybrids, and limiters are proposed as **real code** registered under `src/methods/`, not as hidden parameter files.

3. **Quantitative scoring**  
   Verification emits versioned JSON reports. Absolute and relative scores cover order preservation, excess dissipation, shock quality (thickness / overshoot), and robustness (divergence, NaNs, positivity). An invent loop classifies candidates against the classical `persson_av` baseline.

4. **Human- and agent-driven invention loop**  
   Methods are logged in a laboratory notebook ([`research/experiment_log.md`](research/experiment_log.md)), summarized by frontier/lessons analytics, and—when short-listed—frozen into reproducibility snapshots. The invent scheme (**Gauss–Legendre + Rusanov + SSP-RK3**) is frozen so composite scores remain comparable over time unless a re-baseline is explicitly logged.

Implementation blueprints, residual-hook contracts, and JSON schemas live in [`docs/design.md`](docs/design.md). Publication-oriented VTU/ParaView figures (Schlieren, pressure, sensor, line-outs) are documented in [`docs/visualization.md`](docs/visualization.md); representative Riemann cfg6 and reduced Double-Mach images appear under [Test cases](#test-cases) below (local documentation path; not required CI).

---

## Numerical schemes

### Spatial discretization

FRForge uses **Flux Reconstruction** on tensor-product elements (1D intervals; 2D quads, including isoparametric curved meshes).

| Axis | Default (invent / scoring) | Also available |
|------|----------------------------|----------------|
| Solution points | **Gauss–Legendre (GL)** | Gauss–Lobatto–Legendre (GLL) |
| Correction | Huynh \(g_{DG}\) (Legendre/Radau construction) | — (core fixed) |
| Interface flux | **Rusanov** (local Lax–Friedrichs) | **HLLC** (Euler) |
| Time | **SSP-RK3** | SSP-RK2 |

On GL nodes with the \(g_{DG}\) correction functions, the scheme recovers a DG-equivalent FR formulation in the sense of Huynh. GLL points are supported for scheme-axis studies and robustness exploration; invent composite history stays on the frozen GL + Rusanov + SSP-RK3 defaults.

**References (spatial / FR):**

- H. T. Huynh, *A flux reconstruction approach to high-order schemes including discontinuous Galerkin methods*, AIAA 2007-4079.
- Related FR/CPR literature: Wang & Gao; Vincent, Castonguay & Jameson (energy-stable FR), among others.

### Interface flux and time integration

- **Rusanov** (local Lax–Friedrichs) is the default numerical flux for invent comparisons—simple, robust, and standard in FR/DG shock-capturing studies.
- **HLLC** (Toro) is available for Euler when exploring less dissipative interfaces; it is not the frozen invent default.
- **SSP-RK3** (Gottlieb–Shu) is the default time integrator; SSP-RK2 is available. Order studies use fixed \(\Delta t\) or carefully controlled CFL so temporal error does not pollute spatial rates.

**References (flux / time):**

- V. V. Rusanov, *Calculation of interaction of non-steady shock waves with obstacles*, J. Comput. Math. Phys. USSR, 1961.
- E. F. Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics*, Springer.
- S. Gottlieb & C.-W. Shu, *Total variation diminishing Runge–Kutta schemes*, Math. Comp. 67 (1998); S. Gottlieb, C.-W. Shu & E. Tadmor, *Strong stability-preserving high-order time discretization methods*, SIAM Rev. 43 (2001).

### Reference capturing method

The **honest baseline** for invention is classical **Persson-style modal sensor + element-local artificial viscosity** (`persson_av`). It is deliberately **mild and classical**—a transparent reference, not a heavily tuned champion—so invent composite margins are easy to interpret:

- A modal (Legendre) smoothness indicator \(\sigma\) per element.
- Element viscosity \(\varepsilon \propto c_{av}\,\sigma\,(h/p)\,\lambda_{\max}\), applied through a BR0-style or local second-difference dissipation path (1D and 2D tensor-product implementations).

An example inventable variant (`scaled_persson`) ships as a structural composition with different defaults—not as the research goal. `NullCapturing` (no sensor/AV) is used for smooth order studies and for excess-dissipation comparisons against an undissipated reference on the same mesh.

**References (capturing):**

- P.-O. Persson & J. Peraire, *Sub-cell shock capturing for discontinuous Galerkin methods*, AIAA 2006-112.
- Related artificial-viscosity / sensor work in DG/FR (e.g. Barter & Darmofal; later modal-sensor refinements).

The invent interface is intentionally wider than this baseline: new sensors, non-element-local operators, hybrid FR–WENO-style traces, and limiters can all be expressed through the hook pipeline without forking the FR residual.

---

## Test cases

Evaluation stresses different failure modes. **Coarse/CI-light meshes** are used for fast screening and required CI; short-listed methods are expected to be confirmed on **finer presentation meshes** (and optional full/nightly suites) before any strong claim.

### Smooth order and conservation

| Case | What it stresses |
|------|------------------|
| 1D linear advection (periodic) | Formal spatial order; mass conservation |
| 1D Euler density wave | Smooth Euler order with capturing off; multi-component conservation |
| 2D Euler density wave / advection | Multi-D tensor-product order |
| Isentropic vortex (Yee-type) | Smooth multi-D Euler with exact solution |
| Periodic mass/momentum/energy | Conservation residuals near machine precision (periodic) |

### 1D discontinuous

| Case | What it stresses |
|------|------------------|
| **Sod shock tube** | Shock, contact, rarefaction; overshoot; shock thickness; positivity |
| **Shu–Osher** shock–entropy wave | High-frequency post-shock entropy waves—classic excess-dissipation test |

Sod metrics can be compared to the exact Riemann solution; Shu–Osher uses a self-converged fine-grid reference (hashed CSV under `test/data/`) rather than an external black-box table alone. With the present AV + explicit SSP-RK settings, Shu–Osher is **most stable at lower \(p\)** (the quant suite uses \(p=1\) as the robust high-order start); higher-\(p\) Shu–Osher remains a **stress test**, not a routine CI gate.

### 2D discontinuous

| Case | What it stresses |
|------|------------------|
| **Lax–Liu 2D Riemann** (esp. **cfg 6**) | Multi-wave interactions, shears/contacts on a Cartesian mesh |
| **Reduced Double-Mach-like** | Inclined shock + reflecting wall (optional/full tier; reduced strength for CI cost) |
| **Reduced forward-facing step** | Solid wall / step geometry (optional tier) |
| Element jumps on curved quads | Capturing + geometry together |

**Lax–Liu configuration 6** (unit square, multi-wave contacts/shears at uniform pressure) is the primary multi-D discontinuous gate. The figures below are a documentation-mesh baseline with classical Persson AV (`persson_av`; presentation mesh \(192\times 192\), \(p=2\), \(t=0.15\))—not invent/CI-light resolution.

<p align="center">
  <img src="docs/images/riemann_cfg6_schlieren.png" alt="Riemann cfg6 numerical Schlieren" width="48%" />
  <img src="docs/images/riemann_cfg6_pressure.png" alt="Riemann cfg6 pressure" width="48%" />
</p>

*Figure: 2D Riemann problem (Lax–Liu cfg 6), FR \(p=2\) with Persson AV. **Left:** numerical Schlieren \(\lvert\nabla\rho\rvert\) (white→black). **Right:** pressure field. Generated via high-order VTU + ParaView tessellate/resample ([`docs/visualization.md`](docs/visualization.md)).*

**Reduced Double-Mach-like reflection** (inclined shock + reflecting wall; reduced strength for HO stability on moderate meshes) exercises solid/reflecting BCs and oblique shocks. Baseline below: `persson_av`, presentation mesh \(280\times 100\), \(p=1\), \(t=0.08\).

<p align="center">
  <img src="docs/images/double_mach_schlieren.png" alt="Double Mach numerical Schlieren" width="90%" />
</p>

<p align="center">
  <img src="docs/images/double_mach_pressure.png" alt="Double Mach pressure" width="90%" />
</p>

*Figure: Reduced Double-Mach-like configuration. **Top:** numerical Schlieren. **Bottom:** pressure. Same documentation pipeline as Riemann; optional/full CI tier for the short reduced case.*

### Supporting checks

| Check | Role |
|-------|------|
| Positivity of density/pressure | Hard gate on Euler discontinuous runs |
| Freestream preservation | Cartesian and **curved isoparametric** quads (metric GCL-style residual) |
| BC freestream (transmissive / Dirichlet) | Non-periodic boundary path |
| Robustness matrix | Scheme/method combinations beyond the invent default (`frforge robustness`) |

**Classic references (problems):**

- G. A. Sod, *A survey of several finite difference methods for systems of nonlinear hyperbolic conservation laws*, J. Comput. Phys. 27 (1978).
- C.-W. Shu & S. Osher, *Efficient implementation of essentially non-oscillatory shock-capturing schemes, II*, J. Comput. Phys. 83 (1989).
- P. D. Lax & X.-D. Liu, *Solution of two-dimensional Riemann problems of gas dynamics by positive schemes*, SIAM J. Sci. Comput. 19 (1998).
- P. Woodward & P. Colella, *The numerical simulation of two-dimensional fluid flow with strong shocks*, J. Comput. Phys. 54 (1984) (Double Mach / step-type configurations).
- H. C. Yee, N. D. Sandham & M. J. Djomehri, *Low-dissipative high-order shock-capturing methods using characteristic-based filters*, J. Comput. Phys. 150 (1999) (isentropic vortex setup commonly used as smooth Euler test).

---

## How methods are judged

### Score components

Each verification report carries scores (formula version 1; weights in the report JSON):

| Score | Weight (default) | Meaning |
|-------|------------------|---------|
| **Order preservation** | 0.30 | Smooth-order cases must pass formal-order gates (binary in absolute scoring) |
| **Dissipation** | 0.25 | Low **excess dissipation** vs exact or `NullCapturing` reference in smooth post-shock regions |
| **Shock quality** | 0.25 | Shock thickness (in solution-point spacings) and overshoot |
| **Robustness** | 0.20 | No divergence / NaNs; positivity where required |
| **Composite** | — | Weighted sum of the four |

Invent comparisons also use **relative** maps vs the baseline report (e.g. less dissipation than `persson_av` scores higher) and an **order-vs-dissipation trade-off** check: a method should not regress order while claiming a win.

### Candidate status

| Status | Meaning |
|--------|---------|
| `rejected` | Hard gates failed (divergence, NaN, order/positivity/conservation failures, etc.) |
| `pass_gates` | Valid run, but composite margin vs baseline is below threshold \(\delta\) (default 0.02) or trade-off failed |
| `promising` | Composite beats baseline by \(\geq \delta\) **and** trade-off OK |
| `accepted_candidate` | `promising` **and** high-order VTK produced for inspection |

Narrative promotion in the experiment log (e.g. toward `publication_grade`) additionally expects written **hypothesis** and **lessons**, robustness evidence, and **fine-mesh confirmation**—not only a green coarse invent JSON.

### Performance / threading

**Threaded residuals are for local documentation runs only; invent composite scores and promotion decisions always use the serial residual.**

Opt-in multi-threaded 2-D residual kernels (face traces, volume residual, modal sensor) via:

```bash
# Make Julia workers available, then enable FRForge residual threads
julia -t 8 --project=. scripts/docs/run_vtu_cases.jl \
  --preset presentation --threads 8 --tag baseline
# or: FRFORGE_THREADS=8 julia -t 8 --project=. scripts/docs/...
```

Defaults stay **single-threaded** (bit-deterministic Phase-4 residual). Threaded runs may differ by floating-point roundoff and must not rewrite invent score history. When residual threads > 1, OpenBLAS is pinned to 1 thread to avoid oversubscription.

### Fine-mesh confirmation

Coarse invent meshes are intentionally CI-friendly: they support rapid iteration and stable composite-score history, but they can crown false winners. After a method reaches `promising` / `accepted_candidate` (or is manually short-listed):

1. **Short-list on coarse invent** — `frforge invent` (unchanged; keeps score history).
2. **Confirm on finer meshes** — `frforge confirm --method <name>` re-runs 2D Riemann cfg 6, reduced Double-Mach, and (by default) isentropic vortex order on a noticeably denser mesh; baseline must also finish on the same mesh.
3. **Only then** consider `publication_grade` or a paper-facing freeze: use `frforge snapshot freeze … --require-confirm`.

Confirm statuses are separate from invent `candidate_status`:

| Confirm status | Meaning |
|----------------|---------|
| `confirmed` | Fine-mesh hard gates + competitiveness vs baseline |
| `confirmation_failed` | Divergence/NaN, order regression, multi-D fail, or baseline did not finish |

Default preset `confirm` targets a ~10–30 min class (e.g. Riemann \(64^2\), \(p=2\); DMR \(120\times 40\)). Use `--preset presentation` for paper-grade meshes (hours) or `--preset quick` for smoke. Invent composite scores are **not** rewritten by confirm.

### Laboratory memory

| Tool | Role |
|------|------|
| [`research/experiment_log.md`](research/experiment_log.md) | Authoritative notebook: hypotheses, metrics, strengths/weaknesses, lessons |
| `frforge log summary \| frontier \| lessons \| show` | Structured analytics (frontier includes `confirmed`; summary cautions `confirmation_failed`) |
| `frforge confirm` | Fine-mesh multi-D gate after short-list |
| `frforge snapshot freeze \| verify \| tables` | Reproducibility packages; prefer `--require-confirm` for papers |

---

## Getting started

### Requirements

- Julia **≥ 1.10** (developed/tested on 1.11.x)
- Git

**CI:** Required GitHub Actions runs on **Ubuntu** (Julia 1.10 and 1.11) for PRs to `develop` / `main`. Required CI uses the default invent scheme, modest grids, and reduced suites (~10–15 min budget). Heavy 2D presentation meshes, full robustness matrices, and large VTU trees are **not** required CI. macOS is a first-class local development platform.

### Setup

```bash
git clone https://github.com/RcktMan77/FRForge.git
cd FRForge
julia --project=. -e 'using Pkg; Pkg.instantiate()'
chmod +x bin/frforge
```

### Essential commands

```bash
# Quantitative suite (order + Sod + Shu–Osher) with scored JSON
./bin/frforge test --suite quant --method persson_av --report results/m5/report.json

# Invent / score a registered method vs Persson baseline (coarse; fast)
./bin/frforge invent --method scaled_persson --baseline persson_av --report-dir results/invent

# After short-list: fine-mesh confirmation (local; not required CI)
./bin/frforge confirm --method scaled_persson --baseline persson_av
# Optional: --preset presentation | --preset quick | --vtk | --no-smooth

# 2D / curved / benchmark light suites
./bin/frforge test --suite 2d --report results/m8/report.json
./bin/frforge test --suite curved --report results/curved/report.json
./bin/frforge test --suite benchmarks --report results/bench/report.json

# Experiment log analytics
./bin/frforge log summary
./bin/frforge log frontier
./bin/frforge log lessons

# High-order VTU for ParaView (≥ 5.5 recommended for Lagrange cells)
./bin/frforge run --case euler_density_wave --p 3 --ne 16 --output results/euler.vtu
```

Full CLI: `./bin/frforge --help`. Unit tests: `julia --project=. -e 'using Pkg; Pkg.test()'`.

### Adding a new capturing method

1. Implement `AbstractCapturingMethod` hooks in `src/methods/MyMethod.jl` (override only what you need).
2. `include` and `register_method!("my_method", …)` in [`src/methods/Registry.jl`](src/methods/Registry.jl).
3. Run `./bin/frforge invent --method my_method --baseline persson_av --report-dir results/invent`.
4. Read `candidate_status` and append a thoughtful entry to the experiment log if the result is interesting (`hypothesis` / `lessons` required for promising+ narrative status).
5. If short-listed: `./bin/frforge confirm --method my_method` (fine mesh), then robustness as needed.
6. Paper-facing freeze: `frforge snapshot freeze … --require-confirm` (warns without confirm; hard-fails with the flag).

Structural novelty (new sensors, operators, hybrids) is the primary research goal; coefficient-only search is allowed but secondary.

### Further documentation

| Document | Content |
|----------|---------|
| [`docs/design.md`](docs/design.md) | Architecture, FR conventions, hooks, JSON schema, scoring contracts, historical milestone plan |
| [`docs/visualization.md`](docs/visualization.md) | Local high-order VTU + ParaView publication figures (Riemann cfg6, reduced DMR) |
| [`docs/reproduce_method_template.md`](docs/reproduce_method_template.md) | Template for reproducing a short-listed method |
| [`research/experiment_log.md`](research/experiment_log.md) | Laboratory notebook and invent history |

### Current capability (compact)

| Area | Status |
|------|--------|
| 1D/2D FR (GL/GLL), Rusanov/HLLC, SSP-RK2/3 | Implemented |
| `persson_av` baseline + invent registry | Implemented |
| Quant suite + invent classification | Implemented |
| HO VTK (discontinuous Lagrange) | Implemented |
| Curved quads + freestream gate | Implemented |
| 2D Riemann, vortex, optional DMR/FFS | Implemented (tiered CI) |
| Experiment log analytics + snapshots | Implemented |
| Fine-mesh `frforge confirm` | Implemented (local/manual; not required CI) |

Detailed milestone chronology lives in [`docs/design.md`](docs/design.md) and the experiment log—not in this overview.

---

## License

[MIT](LICENSE) © 2026 RcktMan77

If you use FRForge in research, please cite the repository and relevant classical references above for the schemes and test problems you rely on. Contributions of new capturing methods with scored invent reports and clear experiment-log lessons are especially welcome.
