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

From the repository root. **`--preset presentation`** uses meshes finer than CI-light
(same \(p\) and \(t_f\) as the documentation baseline). Not used by invent or required CI.

| Case | Presentation mesh | \(p\) | \(t_f\) | CI-light (for comparison) |
|------|-------------------|------|---------|---------------------------|
| Riemann cfg6 | \(80\times 80\) | 2 | 0.15 | \(16\times 16\), \(p=1\), \(t=0.08\) |
| DMR reduced | \(120\times 40\) | 1 | 0.08 | \(16\times 4\), \(t=0.03\) |

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

The script **tessellates** HO Lagrange cells, computes Schlieren on the refined field,
uses percentile-based contrast, tight framing, and high-resolution PNGs.

```bash
# macOS — ParaView app pvpython (tested: 6.1.0)
export PVPYTHON="/Applications/ParaView-6.1.0.app/Contents/bin/pvpython"

# Riemann
$PVPYTHON scripts/docs/paraview/plot_2d_publication.py \
  --vtu results/docs_vtu/riemann_cfg6_baseline.vtu \
  --outdir results/docs_figures \
  --prefix riemann_cfg6_baseline \
  --case riemann \
  --subdiv 4 --schlieren-pct 99 --res 3

# Double Mach
$PVPYTHON scripts/docs/paraview/plot_2d_publication.py \
  --vtu results/docs_vtu/double_mach_baseline.vtu \
  --outdir results/docs_figures \
  --prefix double_mach_baseline \
  --case dmr \
  --subdiv 4 --schlieren-pct 99 --res 3
```
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
| `riemann_cfg6_baseline_schlieren.png` | Numerical Schlieren \(\lvert\nabla\rho\rvert\), inverted grayscale |
| `riemann_cfg6_baseline_pressure.png` | Pressure field |
| `riemann_cfg6_baseline_sensor.png` | Modal sensor \(\sigma\) (Cell/Point `sensor`) |
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

## Notes

- **CI:** do not run these cases or `pvpython` in required CI.
- **Resolution:** driver defaults are denser than CI-light gates; edit flags for higher quality (cost grows quickly in 2D HO).
- **Color:** Schlieren uses inverted grayscale; pressure/sensor use colorblind-friendly presets when available in your ParaView build (`Cool to Warm`, `Plasma`).
- **Schlieren in pure ParaView GUI:** open VTU → *Gradient* of `rho` → *Calculator* `mag(Gradients)` → grayscale LUT inverted.
