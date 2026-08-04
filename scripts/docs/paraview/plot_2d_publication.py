#!/usr/bin/env python3
"""
Publication-quality 2D figures from FRForge high-order VTU (ParaView / pvpython).

Tested with ParaView 6.1 (macOS). Tessellates Lagrange HO cells, computes Schlieren
on the refined field, tight framing, explicit scalar ranges, high-res PNG.

  /Applications/ParaView-6.1.0.app/Contents/bin/pvpython \\
    scripts/docs/paraview/plot_2d_publication.py \\
    --vtu results/docs_vtu/riemann_cfg6_baseline.vtu \\
    --outdir results/docs_figures \\
    --prefix riemann_cfg6_baseline \\
    --case riemann

Outputs:
  {prefix}_schlieren.png, _pressure.png, _sensor.png, _lineout.png, _lineout.csv
"""

from __future__ import print_function

import argparse
import csv
import math
import os
import sys
from collections import OrderedDict


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="FRForge 2D publication plots from VTU")
    p.add_argument("--vtu", required=True, help="Input .vtu path")
    p.add_argument("--outdir", required=True, help="Directory for PNG/CSV outputs")
    p.add_argument("--prefix", required=True, help="Filename prefix")
    p.add_argument(
        "--case",
        choices=("riemann", "dmr", "auto"),
        default="auto",
        help="Line-out geometry + aspect presets",
    )
    p.add_argument(
        "--subdiv",
        type=int,
        default=4,
        help="Tessellate MaximumNumberofSubdivisions (HO refinement)",
    )
    p.add_argument(
        "--schlieren-pct",
        type=float,
        default=99.0,
        help="Upper percentile for Schlieren color range (contrast)",
    )
    p.add_argument("--width", type=int, default=0, help="Screenshot width (0=auto)")
    p.add_argument("--height", type=int, default=0, help="Screenshot height (0=auto)")
    p.add_argument("--res", type=int, default=3, help="Screenshot supersampling factor")
    p.add_argument(
        "--no-tessellate",
        action="store_true",
        help="Skip HO tessellation (debug)",
    )
    p.add_argument(
        "--colorbar",
        action="store_true",
        help="Show scalar bars (default: off for tight README crops; ranges print to stdout)",
    )
    return p.parse_args(argv)


def _import_paraview():
    try:
        import paraview.simple  # noqa: F401

        return True
    except Exception as e:
        print(
            "ERROR: import paraview.simple failed. Use ParaView's pvpython.\n  %s" % e,
            file=sys.stderr,
        )
        return False


def _point_array_names(src):
    try:
        info = src.GetPointDataInformation()
        return [info.GetArray(i).GetName() for i in range(info.GetNumberOfArrays())]
    except Exception:
        return []


def _cell_array_names(src):
    try:
        info = src.GetCellDataInformation()
        return [info.GetArray(i).GetName() for i in range(info.GetNumberOfArrays())]
    except Exception:
        return []


def _array_range(src, name, assoc="POINTS"):
    """Return (min, max) for a point/cell array, or None."""
    try:
        info = (
            src.GetPointDataInformation()
            if assoc == "POINTS"
            else src.GetCellDataInformation()
        )
        for i in range(info.GetNumberOfArrays()):
            a = info.GetArray(i)
            if a.GetName() == name:
                r = a.GetRange(0)
                return float(r[0]), float(r[1])
    except Exception:
        pass
    return None


def _fetch_point_scalar_samples(src, name, max_samples=250000):
    """Sample point scalars via servermanager.Fetch for percentile estimates."""
    from paraview import servermanager as sm

    try:
        data = sm.Fetch(src)
        pd = data.GetPointData()
        arr = pd.GetArray(name)
        if arr is None:
            return None
        n = data.GetNumberOfPoints()
        step = max(1, n // max_samples)
        vals = []
        for i in range(0, n, step):
            vals.append(float(arr.GetTuple1(i)))
        return vals
    except Exception as e:
        print("  WARNING: could not fetch samples for %s: %s" % (name, e))
        return None


def _percentile(vals, pct):
    if not vals:
        return None
    s = sorted(vals)
    if pct <= 0:
        return s[0]
    if pct >= 100:
        return s[-1]
    k = (len(s) - 1) * (pct / 100.0)
    f = int(math.floor(k))
    c = min(f + 1, len(s) - 1)
    if f == c:
        return s[f]
    return s[f] + (s[c] - s[f]) * (k - f)


def _setup_view(size):
    from paraview.simple import GetActiveViewOrCreate

    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = list(size)
    view.Background = [1.0, 1.0, 1.0]
    for attr, val in (
        ("OrientationAxesVisibility", 0),
        ("CameraParallelProjection", 1),
        ("UseColorPaletteForBackground", 0),
    ):
        try:
            setattr(view, attr, val)
        except Exception:
            pass
    return view


def _hide_all(view):
    from paraview.simple import GetSources, Hide

    for _, src in GetSources().items():
        try:
            Hide(src, view)
        except Exception:
            pass
    try:
        from paraview.simple import GetScalarBars

        for sb in list(GetScalarBars(view=view)):
            try:
                sb.Visibility = 0
            except Exception:
                pass
    except Exception:
        pass


def _apply_lut(name, presets, invert=False, rng=None, log=False, schlieren_white_to_black=False):
    from paraview.simple import GetColorTransferFunction, GetOpacityTransferFunction

    lut = GetColorTransferFunction(name)
    if schlieren_white_to_black and rng is not None:
        # Explicit: 0 → white, max → black (classic numerical Schlieren)
        lo, hi = float(rng[0]), float(rng[1])
        if hi <= lo:
            hi = lo + 1.0
        try:
            lut.RGBPoints = [
                lo, 1.0, 1.0, 1.0,
                hi, 0.0, 0.0, 0.0,
            ]
            lut.ColorSpace = "RGB"
            try:
                lut.UseLogScale = 0
            except Exception:
                pass
        except Exception:
            for preset in presets:
                try:
                    lut.ApplyPreset(preset, True)
                    lut.InvertTransferFunction()
                    break
                except Exception:
                    continue
            try:
                lut.RescaleTransferFunction(lo, hi)
            except Exception:
                pass
    else:
        for preset in presets:
            try:
                lut.ApplyPreset(preset, True)
                break
            except Exception:
                continue
        if invert:
            try:
                lut.InvertTransferFunction()
            except Exception:
                pass
        if rng is not None:
            lo, hi = float(rng[0]), float(rng[1])
            if hi > lo:
                try:
                    lut.RescaleTransferFunction(lo, hi)
                except Exception:
                    pass
    if rng is not None:
        lo, hi = float(rng[0]), float(rng[1])
        if hi > lo:
            try:
                pwf = GetOpacityTransferFunction(name)
                pwf.RescaleTransferFunction(lo, hi)
            except Exception:
                pass
    if log:
        try:
            lut.MapControlPointsToLogSpace()
            lut.UseLogScale = 1
        except Exception:
            pass
    return lut


def _style_display(
    disp,
    view,
    color_by,
    presets,
    invert=False,
    rng=None,
    show_bar=True,
    schlieren_white_to_black=False,
):
    from paraview.simple import ColorBy

    ColorBy(disp, color_by)
    assoc, name = color_by
    if rng is None:
        try:
            disp.RescaleTransferFunctionToDataRange(True, False)
        except Exception:
            pass
    lut = _apply_lut(
        name,
        presets,
        invert=invert,
        rng=rng,
        schlieren_white_to_black=schlieren_white_to_black,
    )
    try:
        disp.SetScalarBarVisibility(view, bool(show_bar))
    except Exception:
        pass
    try:
        from paraview.simple import GetScalarBar

        sb = GetScalarBar(lut, view)
        sb.Title = name
        sb.ComponentTitle = ""
        sb.Visibility = 1 if show_bar else 0
        try:
            sb.ScalarBarLength = 0.5
            sb.ScalarBarThickness = 14
            sb.TitleFontSize = 14
            sb.LabelFontSize = 11
        except Exception:
            pass
    except Exception:
        pass


def _tight_camera(view, bounds, view_size, margin=0.02):
    """Orthographic XY view tightly framed — domain fills the image (minimal whitespace)."""
    xmin, xmax, ymin, ymax = bounds[0], bounds[1], bounds[2], bounds[3]
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    dx = max(xmax - xmin, 1e-12)
    dy = max(ymax - ymin, 1e-12)
    vw, vh = float(view_size[0]), float(view_size[1])
    aspect_view = vw / max(vh, 1.0)
    # CameraParallelScale = half of visible *height* in world units
    # visible_width = 2 * scale * aspect_view
    # Fit data: need 2*scale >= dy and 2*scale*aspect_view >= dx
    scale_h = 0.5 * dy
    scale_w = 0.5 * dx / aspect_view
    scale = max(scale_h, scale_w) * (1.0 + margin)
    try:
        view.CameraParallelProjection = 1
        view.CameraPosition = [cx, cy, max(dx, dy) * 5.0]
        view.CameraFocalPoint = [cx, cy, 0.0]
        view.CameraViewUp = [0.0, 1.0, 0.0]
        view.CameraParallelScale = scale
        view.CenterOfRotation = [cx, cy, 0.0]
    except Exception:
        from paraview.simple import ResetCamera

        ResetCamera(view)


def _view_size_for_case(case, bounds, width, height):
    if width > 0 and height > 0:
        return width, height
    xmin, xmax, ymin, ymax = bounds[0], bounds[1], bounds[2], bounds[3]
    dx = max(xmax - xmin, 1e-12)
    dy = max(ymax - ymin, 1e-12)
    aspect = dx / dy
    if case == "dmr":
        # wide domain
        w = 2000
        h = max(int(w / aspect), 500)
    else:
        # square-ish
        w = 1600
        h = 1600
    return w, h


def _crop_whitespace(path, pad=8, thr=252):
    """Tight-crop near-white margins (keeps a small pad)."""
    try:
        from PIL import Image
    except Exception:
        return
    im = Image.open(path).convert("RGB")
    w, h = im.size
    px = im.load()
    minx, miny, maxx, maxy = w, h, 0, 0
    found = False
    for y in range(h):
        for x in range(w):
            r, g, b = px[x, y]
            if r < thr or g < thr or b < thr:
                found = True
                if x < minx:
                    minx = x
                if y < miny:
                    miny = y
                if x > maxx:
                    maxx = x
                if y > maxy:
                    maxy = y
    if not found:
        return
    minx = max(0, minx - pad)
    miny = max(0, miny - pad)
    maxx = min(w - 1, maxx + pad)
    maxy = min(h - 1, maxy + pad)
    # Avoid cropping a thin colorbar-only strip incorrectly: require reasonable area
    if (maxx - minx) * (maxy - miny) < 0.05 * w * h:
        return
    cropped = im.crop((minx, miny, maxx + 1, maxy + 1))
    cropped.save(path, optimize=True)
    print("  cropped → %dx%d" % cropped.size)


def _save(view, path, width, height, res=3, crop=True):
    from paraview.simple import SaveScreenshot, Render

    Render(view)
    try:
        SaveScreenshot(
            path,
            view,
            ImageResolution=[int(width * res), int(height * res)],
            FontScaling="Do not scale fonts",
            OverrideColorPalette="WhiteBackground",
            TransparentBackground=0,
            CompressionLevel="1",
        )
    except TypeError:
        SaveScreenshot(
            path,
            view,
            ImageResolution=[int(width * res), int(height * res)],
        )
    if crop:
        _crop_whitespace(path)
    print("  wrote", path)


def _bounds(src):
    from paraview.simple import UpdatePipeline

    UpdatePipeline(proxy=src)
    return src.GetDataInformation().GetBounds()


def _lineout_presets(case, bounds):
    xmin, xmax, ymin, ymax = bounds[0], bounds[1], bounds[2], bounds[3]
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    if case == "dmr":
        # wall-parallel cut and through the triple-point region
        y1 = ymin + 0.25 * (ymax - ymin)
        return [
            ("y_0.25H", (xmin, y1, 0.0), (xmax, y1, 0.0)),
            ("x_mid", (cx, ymin, 0.0), (cx, ymax, 0.0)),
        ]
    return [
        ("y_0.5", (xmin, cy, 0.0), (xmax, cy, 0.0)),
        ("x_0.5", (cx, ymin, 0.0), (cx, ymax, 0.0)),
    ]


def _export_lineouts(src, lines, csv_path, png_path):
    from paraview.simple import PlotOverLine, UpdatePipeline, Delete
    from paraview import servermanager as sm

    rows = []
    headers = ["line", "s", "x", "y", "rho", "p"]
    for name, p0, p1 in lines:
        pol = PlotOverLine(Input=src)
        pol.Point1 = list(p0)
        pol.Point2 = list(p1)
        try:
            pol.Resolution = 800
        except Exception:
            pass
        UpdatePipeline(proxy=pol)
        try:
            fetcher = sm.Fetch(pol)
            n = fetcher.GetNumberOfPoints()
            pd = fetcher.GetPointData()
            arr_rho = pd.GetArray("rho")
            arr_p = pd.GetArray("p")
            for i in range(n):
                pt = fetcher.GetPoint(i)
                s = float(i) / max(n - 1, 1)
                rho = arr_rho.GetValue(i) if arr_rho is not None else float("nan")
                pr = arr_p.GetValue(i) if arr_p is not None else float("nan")
                rows.append([name, s, pt[0], pt[1], rho, pr])
        except Exception as e:
            print("  WARNING: line-out fetch failed for %s: %s" % (name, e))
        try:
            Delete(pol)
        except Exception:
            pass

    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(headers)
        w.writerows(rows)
    print("  wrote", csv_path, "(%d samples)" % len(rows))

    if not rows:
        return
    _lineout_png_matplotlib(rows, png_path)


def _lineout_png_matplotlib(rows, png_path):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print("  WARNING: matplotlib unavailable for lineout PNG:", e)
        return False

    series = OrderedDict()
    for name, s, x, y, rho, p in rows:
        series.setdefault(name, {"s": [], "rho": [], "p": []})
        series[name]["s"].append(s)
        series[name]["rho"].append(rho)
        series[name]["p"].append(p)

    n = len(series)
    fig, axes = plt.subplots(n, 1, figsize=(9.5, 3.2 * n), squeeze=False, dpi=160)
    for ax, (name, data) in zip(axes[:, 0], series.items()):
        ax.plot(data["s"], data["rho"], color="k", lw=2.2, label=r"density $\rho$")
        ax.plot(
            data["s"],
            data["p"],
            color="#0072B2",
            lw=2.0,
            ls="--",
            label=r"pressure $p$",
        )
        ax.set_xlabel(r"normalized arc length $s$", fontsize=11)
        ax.set_ylabel(name, fontsize=11)
        ax.legend(frameon=False, loc="best", fontsize=10)
        ax.grid(True, alpha=0.35, lw=0.6)
        ax.set_title("line-out: %s" % name, fontsize=12)
        # tight y from data with small pad
        ys = [v for v in data["rho"] + data["p"] if v == v]
        if ys:
            lo, hi = min(ys), max(ys)
            pad = 0.05 * (hi - lo + 1e-12)
            ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlim(0.0, 1.0)
        ax.tick_params(labelsize=10)
    fig.tight_layout(pad=0.6)
    fig.savefig(png_path, dpi=180, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("  wrote", png_path)
    return True


def _tessellate(reader, subdiv):
    from paraview.simple import Tessellate, UpdatePipeline

    t = Tessellate(Input=reader)
    try:
        t.MaximumNumberofSubdivisions = int(subdiv)
    except Exception:
        pass
    try:
        t.ChordError = 1.0e-4
    except Exception:
        pass
    try:
        t.MergePoints = 1
    except Exception:
        pass
    try:
        t.OutputDimension = 2
    except Exception:
        # 2 may not be valid on all builds; leave default
        pass
    UpdatePipeline(proxy=t)
    return t


def main(argv=None):
    if not _import_paraview():
        return 2

    from paraview.simple import (
        XMLUnstructuredGridReader,
        Gradient,
        Calculator,
        Show,
        UpdatePipeline,
        Delete,
    )

    args = parse_args(argv)
    vtu = os.path.abspath(args.vtu)
    outdir = os.path.abspath(args.outdir)
    prefix = args.prefix
    os.makedirs(outdir, exist_ok=True)

    if not os.path.isfile(vtu):
        print("ERROR: VTU not found:", vtu, file=sys.stderr)
        return 1

    print("Reading", vtu)
    reader = XMLUnstructuredGridReader(FileName=[vtu])
    # Load all arrays
    try:
        reader.CellArrayStatus = reader.CellArrayStatus
        reader.PointArrayStatus = reader.PointArrayStatus
    except Exception:
        pass
    UpdatePipeline(proxy=reader)
    point_names = _point_array_names(reader)
    cell_names = _cell_array_names(reader)
    print("  PointData:", point_names)
    print("  CellData:", cell_names)
    if "rho" not in point_names:
        print("ERROR: PointData 'rho' required", file=sys.stderr)
        return 1

    case = args.case
    if case == "auto":
        base = (os.path.basename(vtu) + " " + prefix).lower()
        case = "dmr" if ("double_mach" in base or "dmr" in base) else "riemann"

    # HO tessellation so Lagrange interiors are visible
    if args.no_tessellate:
        source = reader
        print("Tessellate: skipped")
    else:
        print("Tessellate: subdivisions=%d" % args.subdiv)
        source = _tessellate(reader, args.subdiv)
        print("  refined PointData:", _point_array_names(source))

    bounds = _bounds(source)
    width, height = _view_size_for_case(case, bounds, args.width, args.height)
    print("View size: %dx%d  bounds: %s" % (width, height, bounds))
    view = _setup_view((width, height))

    # Data ranges for pressure (clamp negatives for display — numerical overshoots)
    prange = _array_range(source, "p", "POINTS")
    if prange:
        plo = max(prange[0], 0.0) if prange[0] < 0 else prange[0]
        # Soft-clip high outliers via samples if available
        psamples = _fetch_point_scalar_samples(source, "p")
        if psamples:
            phi = _percentile(psamples, 99.5)
            plo2 = _percentile(psamples, 0.5)
            plo = max(plo, plo2) if plo2 is not None else plo
            if plo2 is not None and plo2 < 0:
                plo = 0.0
            phi = phi if phi and phi > plo else prange[1]
        else:
            phi = prange[1]
        prange = (plo, phi)
        print("  pressure display range:", prange)

    # ---- Pressure ----
    print("Pressure…")
    _hide_all(view)
    disp = Show(source, view)
    disp.Representation = "Surface"
    try:
        disp.InterpolateScalarsBeforeMapping = 1
    except Exception:
        pass
    _style_display(
        disp,
        view,
        ("POINTS", "p"),
        ["Cool to Warm", "Viridis (matplotlib)", "Blue to Red Rainbow"],
        rng=prange,
        show_bar=args.colorbar,
    )
    _tight_camera(view, bounds, (width, height), margin=0.02)
    _save(view, os.path.join(outdir, "%s_pressure.png" % prefix), width, height, args.res)

    # ---- Schlieren |∇ρ| on tessellated field ----
    print("Schlieren…")
    _hide_all(view)
    grad = Gradient(Input=source)
    try:
        grad.ScalarArray = ["POINTS", "rho"]
    except Exception as e:
        print("  WARNING Gradient.ScalarArray:", e)
    try:
        grad.ComputeGradient = 1
    except Exception:
        pass
    try:
        grad.ResultArrayName = "Gradient"
    except Exception:
        pass
    UpdatePipeline(proxy=grad)

    gnames = _point_array_names(grad)
    gvec = None
    for cand in ("Gradient", "Gradients", "rhoGradient", "Result"):
        if cand in gnames:
            gvec = cand
            break
    if gvec is None:
        try:
            info = grad.GetPointDataInformation()
            for i in range(info.GetNumberOfArrays()):
                a = info.GetArray(i)
                if a.GetNumberOfComponents() >= 3:
                    gvec = a.GetName()
                    break
        except Exception:
            pass

    if gvec is None:
        print("  ERROR: no gradient array; names=%s" % gnames, file=sys.stderr)
    else:
        calc = Calculator(Input=grad)
        calc.ResultArrayName = "schlieren"
        # magnitude of density gradient
        calc.Function = "mag(%s)" % gvec
        UpdatePipeline(proxy=calc)

        # Contrast: use [0, Pth percentile] rather than absolute max outliers
        samples = _fetch_point_scalar_samples(calc, "schlieren")
        hi = _percentile(samples, args.schlieren_pct) if samples else None
        lo = 0.0
        if hi is None or not (hi > lo):
            rr = _array_range(calc, "schlieren", "POINTS")
            hi = rr[1] if rr else 1.0
        # Avoid degenerate range
        if hi <= lo:
            hi = lo + 1.0
        print("  schlieren color range: [%.4g, %.4g] (pct=%.1f)" % (lo, hi, args.schlieren_pct))

        _hide_all(view)
        d2 = Show(calc, view)
        d2.Representation = "Surface"
        try:
            d2.InterpolateScalarsBeforeMapping = 1
        except Exception:
            pass
        _style_display(
            d2,
            view,
            ("POINTS", "schlieren"),
            ["Grayscale", "X Ray"],
            rng=(lo, hi),
            show_bar=args.colorbar,
            schlieren_white_to_black=True,
        )
        _tight_camera(view, bounds, (width, height), margin=0.02)
        _save(
            view,
            os.path.join(outdir, "%s_schlieren.png" % prefix),
            width,
            height,
            args.res,
        )
        try:
            Delete(calc)
            Delete(grad)
        except Exception:
            pass

    # ---- Sensor / AV ----
    print("Sensor / AV…")
    field = None
    assoc = "POINTS"
    # Prefer point-expanded sensor from FRForge writer
    if "sensor" in _point_array_names(source):
        field, assoc = "sensor", "POINTS"
    elif "sensor" in _cell_array_names(source):
        field, assoc = "sensor", "CELLS"
    elif "av" in _point_array_names(source):
        field, assoc = "av", "POINTS"
    elif "av" in _cell_array_names(source):
        field, assoc = "av", "CELLS"

    if field is None:
        print("  WARNING: no sensor/av in VTU")
    else:
        srange = _array_range(source, field, assoc)
        # Soft upper clip at 99th percentile for visibility of weak activations
        if assoc == "POINTS":
            samples = _fetch_point_scalar_samples(source, field)
            if samples:
                hi = _percentile(samples, 99.5)
                lo = 0.0
                if hi and hi > lo:
                    srange = (lo, hi)
        print("  %s range:" % field, srange)
        _hide_all(view)
        d3 = Show(source, view)
        d3.Representation = "Surface"
        try:
            d3.InterpolateScalarsBeforeMapping = 1
        except Exception:
            pass
        _style_display(
            d3,
            view,
            (assoc, field),
            ["Plasma (matplotlib)", "Cool to Warm", "Rainbow Desaturated"],
            rng=srange,
            show_bar=args.colorbar,
        )
        _tight_camera(view, bounds, (width, height), margin=0.02)
        _save(
            view,
            os.path.join(outdir, "%s_sensor.png" % prefix),
            width,
            height,
            args.res,
        )

    # ---- Line-outs on tessellated field ----
    print("Line-outs…")
    lines = _lineout_presets(case, bounds)
    _export_lineouts(
        source,
        lines,
        os.path.join(outdir, "%s_lineout.csv" % prefix),
        os.path.join(outdir, "%s_lineout.png" % prefix),
    )

    print("Done →", outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
