# Phase 4 timing summary

Measured on macOS local (same machine, warm Julia session where noted).
Git before: `b54bd75` (post P3.3b). After: this P4 commit.

## Focused suite wall times (local)

| Suite | Before P4 (s) | After P4 (s) | Speedup |
|-------|---------------|--------------|---------|
| p33a benchmarks | 5.437 | 3.123 | **1.74×** |
| p33b optional | 0.422 | 0.393 | 1.07× |
| p32 curved | 2.804 | 1.97 | **1.42×** |
| isentropic vortex order | 2.87 | 0.798 | **3.60×** |
| euler2d smooth order | 74.037 | 20.43 | **3.62×** |

## Residual microbenchmark (8×8, p=2, NullCapturing)

| Metric | Before | After |
|--------|--------|-------|
| Bytes / residual | ~458 KiB | ~70 KiB (~6.5× less) |
| Vortex full run alloc | ~8–9 GiB | ~0.9 GiB |

## Numeric fidelity

- Freestream residual still ≲ 1e-12 (affine Cartesian).
- Isentropic vortex L2 (p=2, n=8, t=0.5): `0.02069688015968434` → `0.020696880159684756` (relative ~2e-14).
- Observed order: `2.73226063813189` → `2.73226063813174` (FP noise).
- Residual is bit-deterministic (two calls → identical `du`).

## Full local `Pkg.test` wall time

| | Wall (s) | Notes |
|--|----------|--------|
| Before P4 (post–P3.2, ~599 tests) | ~492 | Pre–P3.3a/b content; includes compile |
| After P4 (post–P3.3a/b + perf, ~637 tests) | ~334 | More tests, still faster overall |

## CI

Full package test wall times: see GHA logs on PRs #33–#35 (Julia 1.10/1.11 Ubuntu).
Target required CI remains ~10–15 min; P4 should only reduce wall time.
