# Verification reference data

## shu_osher_ref.csv

| Field | Value |
|-------|-------|
| Problem | Shu–Osher on [0,10], interface x=1, t=1.8 |
| Generator | PerssonAV FR (p=1, Ne=200, c_av=0.5, CFL=0.1, SSP-RK3) |
| Columns | `x, rho` |
| SHA-256 | `4715235eab8435da1def046e5a81be8a79f1c7bb828c2e71d6abe1c6579a8524` |
| ρ range | [0.8000087701561813, 4.477552082786808] |

### Policy (design OQ1)

1. **Self-generated** fine-grid capturing run (not external-only tables).
2. **Self-convergence:** re-run generator at Ne=100,200,400; check shock locus and density extrema stabilize (see script).
3. **Sanity vs literature:** primary shock near x≈2.4; post-shock high-frequency density waves; ρ roughly O(1)–O(4).
4. NullCapturing alone is often unstable on Shu–Osher at practical resolutions; frozen reference uses the **trusted Persson AV** baseline.
5. Regenerate: `julia --project=. test/data/generate_shu_osher_reference.jl`
