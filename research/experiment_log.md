# FRForge Experiment Log

**Authority:** This file is the laboratory notebook and **authoritative memory** for agents and humans.  
**Rule:** Always **read this log before proposing a new capturing method**. Append after every invent / robustness evaluation.

**Frozen invent scheme** (composite-score history): **GL + Rusanov + SSP-RK3** (`DEFAULT_SCHEME`).  
Do not change defaults for invent comparisons without a **logged re-baseline** entry.  
Configurable axes (P2.2): GLL, HLLC, SSP-RK2 — for robustness / exploration only, not invent score history.

**CI policy (one line):** Every addition declares **required CI** vs **full/nightly/manual**; required CI stays under ~10–15 min on Ubuntu; no large VTU/invent trees as PR artifacts.

**Phase 3.1:** 2D Persson AV path on Cartesian quads (tensor-product modal sensor + BR0 AV). Suite: `frforge test --suite 2d_capturing` (CI-light).

**Phase 3.2:** Curved isoparametric quads + metric residual; free-stream preservation merge gate. Suite: `frforge test --suite curved`.

**Phase 3.3a:** Core multi-D benchmarks (CI-light): isentropic vortex order + 2D Riemann cfg6. Suite: `frforge test --suite benchmarks`.

**Phase 3.3b:** Optional reduced DMR/FFS (full/nightly): reflecting/Ghost BCs + solid mask. Suite: `frforge test --suite optional2d`.

**Phase 4:** Residual buffer reuse + in-place fluxes; invent scheme still GL+Rusanov+SSP-RK3. Outputs stable to FP noise.

---

## Schema

Each entry uses the following fields (Markdown bullets under an `### id` heading):

| Field | Required | Notes |
|-------|----------|--------|
| `id` | always | `YYYYMMDD-method_name-short` unique |
| `date` | always | ISO-8601 |
| `method` | always | Registry name |
| `baseline` | invent runs | Usually `persson_av` |
| `hypothesis` | **required if status ≥ promising** | Why try this |
| `scheme.points` / `flux` / `time` | always | Defaults: GL / Rusanov / SSP-RK3 |
| `metrics.*` | invent/score | composite, scores, `candidate_status` |
| `strengths` / `weaknesses` | recommended | |
| `lessons` | **required if status ≥ promising** | Highest-value agent content |
| `status` | always | `open` \| `shortlisted` \| `robustness_pending` \| `publication_grade` \| `archived` \| `baseline` |
| `artifacts` | when available | Paths to JSON reports |
| `git_ref` | optional | SHA or branch |

**Promotion / narrative rule:** Methods at `promising` or higher (`accepted_candidate`, `publication_grade`) must have non-empty **`hypothesis`** and **`lessons`** (not invent placeholders) before any publication-grade claim. Phase 2 robustness evidence is also required for `publication_grade`.

**Write policy:** Append-only (except typo fixes). Status changes = new entry or a dated status note under the entry.

Optional machine-readable index: [`experiment_log.yaml`](experiment_log.yaml).

---

## Entries

### 20260803-m0-m8-platform

- **date:** 2026-08-03
- **method:** _(platform milestone, not a capturing method)_
- **baseline:** n/a
- **hypothesis:** Build a green-field 1D→2D FR laboratory with pluggable capturing, JSON scoring, and invent/score loop before method invention at scale.
- **scheme:** points=GL, flux=Rusanov, time=SSP-RK3
- **metrics:** n/a (milestones 0–8 complete on main)
- **strengths:** First-principles FR; hook pipeline; quant suite (Sod, Shu–Osher); invent classification; HO VTK; 2D Cartesian FR.
- **weaknesses:** No cumulative experiment log until Phase 2.1; scheme axes fixed; Persson AV primarily 1D; no curved geometry; no robustness matrix.
- **lessons:** (1) Fixed Δt for order studies so temporal error does not pollute spatial order. (2) BR0-style AV needs conservative \(c_{av}\) (e.g. 0.1) to avoid blowup. (3) Modal sensor needs careful treatment at element interiors vs jumps. (4) 2D Euler order can be pre-asymptotic on coarse grids. (5) Invent scoring must use NullCapturing-relative excess dissipation for fair comparison.
- **status:** archived
- **artifacts:** none (platform)
- **git_ref:** main @ post-M8 (merge PR #22 develop→main)

### 20260803-persson_av-baseline

- **date:** 2026-08-03
- **method:** persson_av
- **baseline:** n/a (this *is* the classical reference)
- **hypothesis:** Classical Persson-style modal sensor + element artificial viscosity is the honest baseline for invent comparisons, not the research goal.
- **scheme:** points=GL, flux=Rusanov, time=SSP-RK3
- **metrics:**
  - candidate_status: baseline
  - composite: ~0.919 (from invent baseline report; see artifacts)
  - order_preservation / robustness: typically high when AV is off for smooth order cases
- **strengths:** Stable on Sod/Shu–Osher with tuned \(c_{av}\); simple modal sensor; registered and CI-tested.
- **weaknesses:** Excess dissipation in smooth regions if sensor triggers; 1D-centric implementation; parameter-sensitive.
- **lessons:** Keep NullCapturing for smooth order in the quant suite; use Persson only as invent baseline for discontinuous metrics. Document parameters in method_params JSON.
- **status:** baseline
- **artifacts:**
  - method_report: results/invent/baseline_persson_av.json
- **git_ref:** milestone/4–6 lineage on main

### 20260803-scaled_persson-invent

- **date:** 2026-08-03
- **method:** scaled_persson
- **baseline:** persson_av
- **hypothesis:** A structural variant of Persson with different default \(c_{av}\) (and related params) can improve shock metrics without destroying order, serving as a smoke-test invent method for the M6 loop.
- **scheme:** points=GL, flux=Rusanov, time=SSP-RK3
- **metrics:**
  - candidate_status: pass_gates (margin below δ_score for promising on recorded invent run)
  - composite: ~0.924
  - baseline_composite: ~0.919
  - composite_margin: ~+0.005
  - order_preservation: 1.0
  - dissipation: ~0.944
  - shock_quality: ~0.750
  - tradeoff_ok: true
- **strengths:** Slight composite improvement; exercises registry + invent path end-to-end.
- **weaknesses:** Not a large structural novelty (scaled Persson defaults); margin often below promising threshold (δ=0.02); no robustness matrix yet.
- **lessons:** Invent loop works (classify → JSON → status). Coefficient-scale variants alone rarely clear the promising bar — prefer structural hook innovations. Always record scheme; future HLLC/GLL matrix may change ranking.
- **status:** open
- **artifacts:**
  - method_report: results/invent/method_scaled_persson.json
  - compare: results/invent/compare_scaled_persson_vs_persson_av.json
  - baseline_report: results/invent/baseline_persson_av.json
- **git_ref:** milestone/6 invention loop on main

---

## Append protocol

1. Run `frforge invent --method <name> --baseline persson_av` (default scheme only for score history).
2. Invent **auto-appends** a stub entry (see `FRForge.append_experiment_entry!`).
3. If `candidate_status` is `promising` or higher, **edit the entry** to complete `hypothesis` and `lessons` before shortlisting.
4. After invent: if `promising` or higher, run `frforge robustness --method <name> --matrix full` (local/nightly) and fill hypothesis/lessons before `publication_grade`.
5. Required CI may run `matrix=ci` light cells only — not the full 8-cell product.

## Promotion rule (P2.3)

A method may be marked **`publication_grade`** only if **all** hold:

1. Default scheme (GL + Rusanov + SSP-RK3): still `promising` / `accepted_candidate` vs baseline.
2. HLLC cells: no divergence; order preserved; cell status OK.
3. GLL cells: no catastrophic failure.
4. `hypothesis` and `lessons` filled (not placeholders).
5. Robustness summary logged with all cell statuses (`results/robustness/<method>/summary.json`).
