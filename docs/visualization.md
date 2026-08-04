# Publication figures from high-order VTU

Local post-processing for documentation / README images. **Not** part of invent or required CI.

Frozen invent scheme remains **GL + Rusanov + SSP-RK3** (not modified here).

---

## What is in the VTU?

| Field | Always? | Notes |
|-------|---------|--------|
| `rho`, `u`, `v`, `p`, `rho_u`, `rho_v`, `E` | Yes (Euler 2D) | Primitives + conserved at HO Lagrange nodes |
| `sensor` | **Only with diagnostics** | Per-element Persson modal sensor \(\sigma\), expanded to points + CellData |
| `av` | **Only with diagnostics** | Element artificial viscosity \(\varepsilon = c_{av}\sigma(h/p)\lambda_{\max}\) |

Default `write_vtu_high_order` does **not** write sensor/AV (no invent/CI overhead).

Documentation path:

```julia
write_vtu_high_order_with_capturing(path, state, eq, method)
# or CLI:  frforge run … --output out.vtu --vtk-diagnostics
```

Sensor/AV come from the same hooks as the residual (`sense!` + `element_viscosities_2d`). For `NullCapturing`, both are zero.

---

## Cases

1. **2D Riemann (Lax–Liu cfg 6)** — multi-wave contacts/shears on the unit square.  
2. **Reduced Double-Mach-like** (`strength=:reduced`) — inclined shock + reflecting wall.

Recommended baseline for figures: `persson_av` (classical reference).  
For a candidate: same scripts, change `--method` / VTU tag (e.g. `scaled_persson`).

---

## Exact commands (local)

### 1. Generate presentation VTU (Julia)

From the repository root. **`--preset presentation`** uses meshes much finer than CI-light
(same \(p\) and \(t_f\) as the documentation baseline). Not used by invent or required CI.

| Case | Presentation mesh | \(p\) | \(t_f\) | CI-light (for comparison) |
|------|-------------------|------|---------|---------------------------|
| Riemann cfg6 | \(192\times 192\) | 2 | 0.15 | \(16\times 16\), \(p=1\), \(t=0.08\) |
| DMR reduced | \(280\times 100\) | 1 | 0.08 | \(16\times 4\), \(t=0.03\) |

Wall-clock (Apple Silicon, single process, rough): Riemann \(192^2\) can be **~60–90+ min** serial; DMR adds more.
Use `--preset quick` for a fast smoke of the pipeline.

**Threaded residuals are for local documentation runs only; invent composite scores and promotion decisions always use the serial residual.**

```bash
# Optional: multi-threaded residual (needs julia -t N; not invent-score-safe)
julia -t 8 --project=. scripts/docs/run_vtu_cases.jl \
  --preset presentation --threads 8 --tag baseline
```

Long runs print progress every N SSP-RK3 steps (`--progress-every`, auto-enabled on large meshes).

**Why so fine?** FR/DG fields are discontinuous across faces. On a coarse Cartesian mesh,
raw VTU rendering (and even mild HO tessellation) shows horizontal/vertical bands —
especially Schlieren \(\lvert\nabla\rho\rvert\) when gradients hit face jumps. Dense
presentation meshes + ParaView tessellation + **ResampleToImage** make those artifacts
sub-pixel for README figures. This does **not** change invent defaults or CI gates.

ParaView post-processing defaults (see `plot_2d_publication.py`):

- Tessellate HO Lagrange cells: `--subdiv` default **8**, tight chord/field error, `MergePoints=0`
- **ResampleToImage** before pressure/Schlieren/sensor (auto **1600** Riemann / **2000** DMR on longest axis)
- Gradient of density on the continuous resampled field (avoids mesh-aligned face spikes)

```bash
# Baseline (Persson AV) — presentation quality
julia --project=. scripts/docs/run_vtu_cases.jl \
  --preset presentation \
  --method persson_av \
  --outdir results/docs_vtu \
  --tag baseline

# Optional: candidate method (same meshes)
julia --project=. scripts/docs/run_vtu_cases.jl \
  --preset presentation \
  --method scaled_persson \
  --outdir results/docs_vtu \
  --tag scaled_persson

# Faster smoke (coarser): --preset quick
```

VTU includes `rho`, `p`, … plus **`sensor`** and **`av`** via
`write_vtu_high_order_with_capturing` (documentation path only).

| File | Case |
|------|------|
| `results/docs_vtu/riemann_cfg6_baseline.vtu` | Riemann cfg6 |
| `results/docs_vtu/double_mach_baseline.vtu` | Reduced DMR |

### 2. Generate images (`pvpython`)

The script **tessellates** HO Lagrange cells, **resamples** to a fine image grid,
computes Schlieren on the continuous density field, uses percentile-based contrast,
tight framing, and high-resolution PNGs.

```bash
# macOS — ParaView app pvpython (tested: 6.1.0)
export PVPYTHON="/Applications/ParaView-6.1.0.app/Contents/bin/pvpython"

# Riemann (defaults: --subdiv 8, auto image-res 1600)
$PVPYTHON scripts/docs/paraview/plot_2d_publication.py \
  --vtu results/docs_vtu/riemann_cfg6_baseline.vtu \
  --outdir results/docs_figures \
  --prefix riemann_cfg6_baseline \
  --case riemann \
  --schlieren-pct 98 --res 3

# Double Mach (auto image-res 2000)
$PVPYTHON scripts/docs/paraview/plot_2d_publication.py \
  --vtu results/docs_vtu/double_mach_baseline.vtu \
  --outdir results/docs_figures \
  --prefix double_mach_baseline \
  --case dmr \
  --schlieren-pct 98 --res 3
```

If residual banding remains on an older coarse VTU, raise sampling further, e.g.
`--subdiv 10 --image-res 2048` (memory and time grow quickly).

### 3. Candidate method (same pipeline)

```bash
julia --project=. scripts/docs/run_vtu_cases.jl --method scaled_persson --tag scaled_persson
$PVPYTHON scripts/docs/paraview/plot_2d_publication.py \
  --vtu results/docs_vtu/riemann_cfg6_scaled_persson.vtu \
  --outdir results/docs_figures \
  --prefix riemann_cfg6_scaled_persson \
  --case riemann
```

---

## Output filenames (README assets)

### Riemann cfg6 (baseline)

| File | Content |
|------|---------|
| `riemann_cfg6_baseline_schlieren.png` | Numerical Schlieren \(\lvert\nabla\rho\rvert\), white→black |
| `riemann_cfg6_baseline_pressure.png` | Pressure field |
| `riemann_cfg6_baseline_sensor.png` | Modal sensor \(\sigma\) (resampled) |
| `riemann_cfg6_baseline_lineout.png` | Density/pressure cuts (\(y=1/2\), \(x=1/2\)) |
| `riemann_cfg6_baseline_lineout.csv` | Line-out samples |

### Double Mach reduced (baseline)

| File | Content |
|------|---------|
| `double_mach_baseline_schlieren.png` | Schlieren |
| `double_mach_baseline_pressure.png` | Pressure |
| `double_mach_baseline_sensor.png` | Sensor \(\sigma\) |
| `double_mach_baseline_lineout.png` | Cuts (mid-\(y\), mid-\(x\)) |
| `double_mach_baseline_lineout.csv` | Line-out samples |

All written under `--outdir` (e.g. `results/docs_figures/`). Copy into `docs/images/` if you want them tracked for the README (optional; large binaries often stay untracked).

---

## Artifact notes (banding / blockiness)

| Cause | Mitigation in this pipeline |
|-------|----------------------------|
| DG face discontinuities | Tessellate (`MergePoints=0`) + **ResampleToImage** before display/gradient |
| Schlieren on discontinuous \(\rho\) | Gradient after continuous resample |
| Cartesian mesh alignment (Riemann) | Dense presentation mesh (\(192^2\)) so faces are fine vs features |
| Element-wise sensor | Resample sensor to image for README crops |
| Under-sampled HO interior | `--subdiv 8` + tight chord/field error |

Residual weak banding near strong axis-aligned waves can remain (true under-resolved
structure + mild AV oscillations); further mesh refinement is the physics-side fix,
not invent-scheme changes.

---

## Notes

- **CI:** do not run these cases or `pvpython` in required CI.
- **Resolution:** driver defaults are denser than CI-light gates; edit flags for higher quality (cost grows quickly in 2D HO).
- **Color:** Schlieren uses white→black grayscale; pressure/sensor use colorblind-friendly presets when available (`Cool to Warm`, `Plasma`).
- **Schlieren in pure ParaView GUI:** open VTU → *Tessellate* → *Resample to Image* → *Gradient* of `rho` → *Calculator* `mag(Gradient)` → grayscale LUT white→black.
