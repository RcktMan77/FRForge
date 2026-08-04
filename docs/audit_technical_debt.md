# FRForge technical debt, consistency, and safe-performance audit

| Field | Value |
|-------|--------|
| **Date** | 2026-08-04 |
| **Scope** | Audit + cleanup implementation |
| **Out of scope** | Numerics, physics, invent scoring, frozen scheme (GL+Rusanov+SSP-RK3), serial residual bit-determinism, invent redesign |
| **Cleanup status** | Implemented 2026-08-04 (see § Implementation status below) |

---

## Executive summary

FRForge is a coherent green-field lab: module layout matches the design blueprint, Phase-4 residual discipline is real, and invent/confirm/threading policies are documented. The main debt is **scale and consolidation** (very large case/CLI files), **stale narrative** in `docs/design.md`, a few **dead helpers**, **no formatter config**, and **residual-per-call allocations** that are safe-to-attack only with careful buffer reuse (not free wins without tests).

**Must-fix** items are small and behavior-preserving. **Should-fix** includes formatter + low-risk alloc work. **Optional** covers larger refactors.

---

## A. Consistency & formatting

| ID | Severity | Finding | Location | Rationale |
|----|----------|---------|----------|-----------|
| A1 | Must-fix | No JuliaFormatter / EditorConfig in repo | project root | Style drift will grow; 4-space indent dominates but a few files have 2-space outliers |
| A2 | Must-fix | Indent outliers (2-space lines mixed into 4-space codebase) | e.g. `PerssonAV.jl`, `PerssonAV2D.jl`, `Correction.jl`, `SSP_RK3.jl`, `ExactSod.jl`, `Rusanov.jl` (small counts) | Cosmetic; safe reindent with formatter |
| A3 | Should-fix | Public export surface is very large (~90+ exports) with little public/internal split | `src/FRForge.jl` | Maintainability; many symbols are only for tests/CLI |
| A4 | Should-fix | File headers still say “Milestone N” only | `Cases.jl` L1, various invent modules | Comments lag multi-phase reality |
| A5 | Must-fix | `docs/design.md` “Current state” still describes empty repo + first commit only | `docs/design.md` ~Background | Contradicts shipped M0–M8 / P2–P5 code |
| A6 | Should-fix | Residual2D docstring says Phase-4 reuses buffers; call still `zeros`/`similar` for traces/fhat/`u_work`/`σ` every residual | `Residual2D.jl` residual! | Doc slightly oversells; Phase-4 cut *micro*-allocs (fluxes), not face arrays |
| A7 | Optional | `cli/main.jl` (~1230 LOC) holds all subcommands | `src/cli/main.jl` | Split by command would help reviews; pure move risk if done carefully |

---

## B. Stale & dead code

| ID | Severity | Finding | Location | Rationale |
|----|----------|---------|----------|-----------|
| B1 | Must-fix | Allocating `face_trace(...)` wrapper unused; only `face_trace!` used | `Residual2D.jl` L10–16 | Dead API; remove or mark internal + test if kept |
| B2 | Must-fix | `contravariant_fluxes` defined but never called (volume path inlines same math) | `Residual2D.jl` L65–81 | Dead code; safe delete |
| B3 | Should-fix | `l2_error_all` exported but unused outside definition | `SolutionState.jl`, export in `FRForge.jl` | Dead export; keep if planned API, else unexport/remove |
| B4 | Should-fix | `riemann2d_cfg3_ic` / `RIEMANN2D_CFG3` only referenced in Cases2D setup, not in required CI suites | `Cases2D.jl` | Not dead (research path); document as full/nightly only |
| B5 | Must-fix | Design doc stale “no solver code / only README” | `docs/design.md` | Misleading for newcomers |
| B6 | Optional | Experiment log still carries long phase checklist at top | `research/experiment_log.md` | Useful history; optional trim once roadmap frozen |
| B7 | Should-fix | Working tree mixes many features uncommitted | git status | Not dead code; process risk—split commits before merge |

**Not dead (verified used):** `g_DG_endpoints` (tests), `write_report_skeleton` / `load_report` (CLI/tests), `viscous_mass_residual_scale*` (tests/cases), `dissipation_operator` (VTK), `lagrange_basis_matrix` (VTK).

---

## C. Maintainability

| ID | Severity | Finding | Location | Rationale |
|----|----------|---------|----------|-----------|
| C1 | Should-fix | Massive duplicated “case dict” construction | `Cases.jl` (~1240 LOC), `Cases2D.jl` (~1360 LOC) | Safe helper for common keys (`name`, `pass`, `diverged`, …) reduces bugs without logic change |
| C2 | Should-fix | CLI arg parsers duplicated (points/flux/time) | `cli/main.jl` | Extract shared scheme-arg helper |
| C3 | Optional | Confirm / invent / robustness each reimplement report paths + log append patterns | `invent/*` | Shared “write trio + append log” helper |
| C4 | Optional | 1D vs 2D residual and AV codepaths largely parallel | `Residual.jl` / `Residual2D.jl`, `PerssonAV.jl` / `PerssonAV2D.jl` | Unifying abstractions is high judgment / risk |
| C5 | Must-fix | Module comments still “Phase 3.1” only where 3.2/3.3/thread exist | e.g. `PerssonAV2D.jl` L1 | Comment accuracy |
| C6 | Should-fix | Threading helpers (`VolumeScratch2D`, pools) not exported (good) but not documented in design | `Threading.jl` | Add short note to design or visualization doc only |
| C7 | Optional | No `src/FRForge.jl` section comments distinguishing public invent API vs internal FR | `FRForge.jl` | Export hygiene |

---

## D. Safe performance opportunities

**Constraint reminder:** serial residual must stay bit-deterministic; invent/confirm defaults serial; no CFL/scheme/scoring changes.

| ID | Severity | Finding | Location | Risk if done carefully |
|----|----------|---------|----------|-------------------------|
| D1 | Should-fix | **Every** 2D residual allocates `u_work`, 4× traces, 4× fhat, `σ` | `Residual2D.jl` residual! | Medium: preallocate on `SolutionState2D` or residual cache; must `fill!` same order; re-run bit-determinism test |
| D2 | Should-fix | 1D residual allocates `u_work`, `fL`/`fR`, `σ` each call | `Residual.jl` | Same as D1 |
| D3 | Should-fix | Persson sensor rebuilds Vandermonde and does `V \ U` allocating `Matrix` per element every residual | `PerssonAV2D.jl` sense! | Medium: cache `V`/`inv(V)` on operators or method; solve in-place; FP may change slightly if inverse formulation differs—**prefer factorize once, same `\` algorithm`** |
| D4 | Optional | BR0 AV allocates many large face arrays per residual | `PerssonAV2D.jl` apply_dissipation_br0_2d! | High value, more touch surface |
| D5 | Must-fix (pure overhead only if proven unused) | Dead `contravariant_fluxes` / allocating `face_trace` | Residual2D | None if deleted |
| D6 | Optional | Interface flux loops still serial under threading T1 | Residual2D | Already planned T2; not this audit’s cleanup |
| D7 | Should-fix | `zeros` for `du` fill already uses `fill!`; good | Residual2D | No change |
| D8 | Optional | Presentation/docs wall time dominated by steps × residual; progress + threads already help | docs driver | Done in WIP; not debt |

**Do not treat as “safe cleanup” without dedicated PR + bit tests:** changing loop order, fusing interface+volume, OpenBLAS defaults for invent, CFL, or AV formula.

---

## Proposed cleanup plan (for after audit approval)

### PR-A — Must-fix (safe, small)

1. **Docs accuracy**
   - Refresh `docs/design.md` “Current state” to match shipped platform (M0–M8, P2–P5, invent/confirm/threading policy one-liners).
   - Fix stale file headers on `Cases.jl` / `PerssonAV2D.jl` (comment only).
2. **Dead code**
   - Remove unused `contravariant_fluxes`.
   - Remove unused allocating `face_trace` **or** keep only if you want a public debug helper (prefer remove).
3. **Hygiene**
   - Ensure all `src/**/*.jl` end with a newline.
   - Add minimal `.JuliaFormatter.toml` (4-space, same as majority) **without** reformatting entire tree in same PR *or* full format as separate commit.

### PR-B — Should-fix (still low risk, separate PR)

1. Shared `case_dict_base(...)` helper for verification reports (behavior-preserving field set).
2. Shared CLI scheme-arg parser.
3. Unexport or document `l2_error_all`.
4. Residual buffer fields on state (1D/2D) + fill! — **with** bit-determinism tests mandatory.

### Deferred (optional / backlog)

- Split `cli/main.jl` and `Cases*.jl` into multiple files.
- Full BR0 buffer pool.
- Sensor in-place modal transform.
- Export surface reduction / public API doc.
- Thread face fluxes (T2 performance plan).

---

## Verification checklist (for cleanup PRs)

| Check | Method |
|-------|--------|
| Tests pass | `julia --project=. -e 'using Pkg; Pkg.test()'` |
| Serial residual bit-determinism | Existing/new test: two `residual!` calls → `du1 == du2` under `with_serial_residual` |
| Invent scores unchanged | Run `frforge invent --method scaled_persson` (or quant) before/after on same commit mesh; compare composites to FP noise only if no invent path touched |
| CI budget | No new heavy suites; formatter CI optional later |

---

## Severity summary counts

| Bucket | Count (approx.) |
|--------|-----------------|
| Must-fix | ~7 (A1 partial, A2, A5, A6 doc, B1, B2, B5, C5, D5) |
| Should-fix | ~10 |
| Optional / backlog | ~8 |

---

## What this audit deliberately did **not** do

- No residual math or AV changes  
- No invent/score/confirm scoring changes  
- No requirement to land threading or confirm (already in WIP)  
- No large performance experiments  

---

## Implementation status (2026-08-04)

| Item | Status |
|------|--------|
| Dead `contravariant_fluxes` / allocating `face_trace` | **Removed** |
| Design “Current state” + file headers / residual docs | **Updated** |
| Docstring cleanliness (stale Milestone/Phase headers, signature mismatches, suite lists, design Non-Goals / operators / package tree) | **Updated** (2026-08-04 follow-up) |
| `.JuliaFormatter.toml` | **Added** (full-tree reformat deferred to avoid huge noise PR) |
| `l2_error_all` unexported | **Done** (still available as `FRForge.l2_error_all`) |
| Residual workspaces 1D/2D (reuse traces/fhat/σ/u_work) | **Done** |
| Cached `ops.V_legendre` for Persson sensor | **Done** |
| BR0 AV buffer pool on residual workspace | **Done** |
| Face-parallel interior interfaces (when residual threads > 1) | **Done** (serial order preserved at thr=1) |
| `case_report_dict` helper | **Added** (`CaseReport.jl`); gradual call-site migration optional |
| CLI `scheme_from_cli_opts` | **Added** (`cli/scheme_args.jl`) |
| Split `Cases*.jl` / full `main.jl` | **Partial** — helpers extracted; full file split deferred (high churn, no behavior win) |
| Full JuliaFormatter pass on all sources | **Deferred** (optional follow-up PR) |
| Unify 1D/2D residual abstractions | **Deferred** (high judgment / risk of accidental numeric drift) |

**Verification:** `Pkg.test()` — **772 tests passed** after cleanup (2026-08-04). Serial residual bit-equality holds under `with_serial_residual`.

## Next step

Remaining optional churn (full formatter run, Cases/CLI file splits, 1D/2D residual unification) can be separate PRs if still desired.

---

## Appendix: size hotspots

| File | ~LOC | Role |
|------|------|------|
| `Cases2D.jl` | 1360 | 2D verification cases |
| `Cases.jl` | 1240 | 1D verification cases |
| `cli/main.jl` | 1230 | All CLI |
| `Confirm.jl` | 710 | Fine-mesh confirm |
| `LogAnalytics.jl` | 610 | Log parser/views |
| `Residual2D.jl` | 450 | 2D residual |

---

## Appendix: positive findings (keep)

- Phase-4 in-place fluxes and residual arithmetic-order comments are intentional and valuable.  
- Invent forces serial residual; confirm-for-promotion serial by default; docs-only threading is gated.  
- JSON schema validation (`validate_report_keys`) is centralized.  
- Capturing hooks are type-stable dispatch without residual naming concrete methods.  
- Test layout mirrors subsystems (`test_*.jl` includes).  
