#!/usr/bin/env python3
"""
Publication-style 2D figures from FRForge high-order VTU (ParaView / pvpython).

Tested with ParaView 6.1 (macOS). Use the `pvpython` that ships with ParaView:

  /Applications/ParaView-6.1.0.app/Contents/bin/pvpython \\
    scripts/docs/paraview/plot_2d_publication.py \\
    --vtu results/docs_vtu/riemann_cfg6_baseline.vtu \\
    --outdir results/docs_figures \\
    --prefix riemann_cfg6_baseline \\
    --case riemann

Produces:
  {prefix}_schlieren.png   — inverted grayscale |∇ρ|
  {prefix}_pressure.png    — pressure field
  {prefix}_sensor.png      — modal sensor σ (or av if sensor missing)
  {prefix}_lineout.png     — density/pressure line cuts (when chart view works)
  {prefix}_lineout.csv     — numeric samples for the cuts
"""

from __future__ import print_function

import argparse
import csv
import os
import sys


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="FRForge 2D publication plots from VTU")
    p.add_argument("--vtu", required=True, help="Input .vtu path")
    p.add_argument("--outdir", required=True, help="Directory for PNG/CSV outputs")
    p.add_argument(
        "--prefix",
        required=True,
        help="Filename prefix, e.g. riemann_cfg6_baseline",
    )
    p.add_argument(
        "--case",
        choices=("riemann", "dmr", "auto"),
        default="auto",
        help="Line-out geometry presets",
    )
    p.add_argument("--width", type=int, default=1400, help="Screenshot width")
    p.add_argument("--height", type=int, default=1200, help="Screenshot height")
    p.add_argument(
        "--res",
        type=int,
        default=2,
        help="Magnification for SaveScreenshot (supersampling)",
    )
    return p.parse_args(argv)


def _import_paraview():
    try:
        import paraview.simple  # noqa: F401
        return True
    except Exception as e:
        print(
            "ERROR: Could not import ParaView Python bindings.\n"
            "  Run with the ParaView-bundled pvpython, e.g.:\n"
            "    /Applications/ParaView-6.1.0.app/Contents/bin/pvpython \\\n"
            "      scripts/docs/paraview/plot_2d_publication.py ...\n"
            "Import error: %s" % e,
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


def _setup_view(size=(1400, 1200)):
    from paraview.simple import GetActiveViewOrCreate

    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = list(size)
    view.Background = [1.0, 1.0, 1.0]
    try:
        view.OrientationAxesVisibility = 0
    except Exception:
        pass
    try:
        view.CameraParallelProjection = 1
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
    # Remove leftover scalar bars from previous plots (ParaView keeps them otherwise)
    try:
        from paraview.simple import GetScalarBars

        for sb in list(GetScalarBars(view=view)):
            try:
                sb.Visibility = 0
            except Exception:
                pass
    except Exception:
        pass
    try:
        # Alternate API
        bars = getattr(view, "Representations", None)
        if bars is not None:
            for rep in list(view.Representations):
                try:
                    if "ScalarBar" in type(rep).__name__ or "LUT" in type(rep).__name__:
                        rep.Visibility = 0
                except Exception:
                    pass
    except Exception:
        pass


def _apply_lut(name, presets, invert=False):
    from paraview.simple import GetColorTransferFunction

    lut = GetColorTransferFunction(name)
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
    return lut


def _save(view, path, width, height, res=2):
    from paraview.simple import SaveScreenshot, Render

    Render(view)
    kwargs = dict(
        ImageResolution=[int(width * res), int(height * res)],
        TransparentBackground=0,
    )
    # Palette keyword differs slightly across versions
    try:
        SaveScreenshot(path, view, OverrideColorPalette="WhiteBackground", **kwargs)
    except TypeError:
        try:
            SaveScreenshot(path, view, **kwargs)
        except Exception as e:
            print("  ERROR SaveScreenshot:", e, file=sys.stderr)
            return False
    print("  wrote", path)
    return True


def _bounds(reader):
    from paraview.simple import UpdatePipeline

    UpdatePipeline(proxy=reader)
    d = reader.GetDataInformation().GetBounds()
    return d  # xmin,xmax,ymin,ymax,zmin,zmax


def _lineout_presets(case, bounds):
    xmin, xmax, ymin, ymax = bounds[0], bounds[1], bounds[2], bounds[3]
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    if case == "dmr":
        return [
            ("y_mid", (xmin, cy, 0.0), (xmax, cy, 0.0)),
            ("x_mid", (cx, ymin, 0.0), (cx, ymax, 0.0)),
        ]
    return [
        ("y_0.5", (xmin, cy, 0.0), (xmax, cy, 0.0)),
        ("x_0.5", (cx, ymin, 0.0), (cx, ymax, 0.0)),
    ]


def _export_lineouts(reader, lines, csv_path, png_path, width, height, res):
    from paraview.simple import (
        PlotOverLine,
        UpdatePipeline,
        CreateView,
        Show,
        SaveScreenshot,
        Render,
        Delete,
    )
    from paraview import servermanager as sm

    rows = []
    headers = ["line", "s", "x", "y", "rho", "p"]
    first_pol = None

    for name, p0, p1 in lines:
        pol = PlotOverLine(Input=reader)
        pol.Point1 = list(p0)
        pol.Point2 = list(p1)
        try:
            pol.Resolution = 500
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

        if first_pol is None:
            first_pol = pol
        else:
            try:
                Delete(pol)
            except Exception:
                pass

    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(headers)
        w.writerows(rows)
    print("  wrote", csv_path, "(%d samples)" % len(rows))

    if first_pol is None or len(rows) == 0:
        print("  WARNING: no line-out samples; skip PNG")
        return

    # Prefer matplotlib if available (cleaner, version-stable)
    if _lineout_png_matplotlib(rows, png_path, width, height):
        return

    # Fallback: ParaView XY chart
    try:
        chart = CreateView("XYChartView")
        chart.ViewSize = [width, height]
        disp = Show(first_pol, chart)
        try:
            disp.SeriesVisibility = ["rho", "p"]
        except Exception:
            pass
        Render(chart)
        SaveScreenshot(
            png_path,
            chart,
            ImageResolution=[width * res, height * res],
        )
        print("  wrote", png_path)
        try:
            Delete(chart)
        except Exception:
            pass
    except Exception as e:
        print("  WARNING: could not save lineout PNG:", e)


def _lineout_png_matplotlib(rows, png_path, width, height):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return False

    # group by line name
    from collections import OrderedDict

    series = OrderedDict()
    for name, s, x, y, rho, p in rows:
        series.setdefault(name, {"s": [], "rho": [], "p": []})
        series[name]["s"].append(s)
        series[name]["rho"].append(rho)
        series[name]["p"].append(p)

    n = len(series)
    fig_w = max(width / 100.0, 8)
    fig_h = max(height / 100.0, 3 * n)
    fig, axes = plt.subplots(n, 1, figsize=(fig_w, fig_h), squeeze=False)
    for ax, (name, data) in zip(axes[:, 0], series.items()):
        ax.plot(data["s"], data["rho"], "k-", lw=1.5, label=r"$\rho$")
        ax.plot(data["s"], data["p"], "C0--", lw=1.2, label=r"$p$")
        ax.set_xlabel("normalized arc length s")
        ax.set_ylabel(name)
        ax.legend(frameon=False, loc="best")
        ax.grid(True, alpha=0.3)
        ax.set_title("line-out: %s" % name)
    fig.tight_layout()
    fig.savefig(png_path, dpi=150, facecolor="white")
    plt.close(fig)
    print("  wrote", png_path, "(matplotlib)")
    return True


def _top_down_camera(view, bounds):
    """Orthographic top-down view of the XY plane."""
    xmin, xmax, ymin, ymax = bounds[0], bounds[1], bounds[2], bounds[3]
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    dx = max(xmax - xmin, 1e-12)
    dy = max(ymax - ymin, 1e-12)
    try:
        view.CameraParallelProjection = 1
        view.CameraPosition = [cx, cy, max(dx, dy) * 3.0]
        view.CameraFocalPoint = [cx, cy, 0.0]
        view.CameraViewUp = [0.0, 1.0, 0.0]
        view.CameraParallelScale = 0.55 * max(dx, dy)
    except Exception:
        from paraview.simple import ResetCamera

        ResetCamera(view)


def main(argv=None):
    if not _import_paraview():
        return 2

    # Lazy imports after path check — only names present in PV 5.x/6.x simple API
    from paraview.simple import (
        XMLUnstructuredGridReader,
        Gradient,
        Calculator,
        Show,
        UpdatePipeline,
        ColorBy,
        Render,
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
    UpdatePipeline(proxy=reader)
    point_names = _point_array_names(reader)
    cell_names = _cell_array_names(reader)
    print("  PointData:", point_names)
    print("  CellData:", cell_names)

    if "rho" not in point_names:
        print("ERROR: PointData 'rho' required for Schlieren/pressure plots.", file=sys.stderr)
        return 1

    case = args.case
    if case == "auto":
        base = os.path.basename(vtu).lower() + " " + prefix.lower()
        case = "dmr" if ("double_mach" in base or "dmr" in base) else "riemann"

    view = _setup_view((args.width, args.height))
    bounds = _bounds(reader)

    # ---- Pressure ----
    print("Pressure field…")
    _hide_all(view)
    disp = Show(reader, view)
    disp.Representation = "Surface"
    if "p" in point_names:
        ColorBy(disp, ("POINTS", "p"))
        disp.RescaleTransferFunctionToDataRange(True, False)
        disp.SetScalarBarVisibility(view, True)
        _apply_lut("p", ["Cool to Warm", "Viridis (matplotlib)", "Blue to Red Rainbow"])
    _top_down_camera(view, bounds)
    _save(
        view,
        os.path.join(outdir, "%s_pressure.png" % prefix),
        args.width,
        args.height,
        args.res,
    )

    # ---- Schlieren |∇ρ| inverted grayscale (Gradient filter, PV 5/6) ----
    print("Numerical Schlieren…")
    _hide_all(view)
    grad = Gradient(Input=reader)
    try:
        grad.ScalarArray = ["POINTS", "rho"]
    except Exception as e:
        print("  WARNING: could not set Gradient.ScalarArray:", e)
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
        # pick first multi-component array that is not an original scalar
        try:
            info = grad.GetPointDataInformation()
            for i in range(info.GetNumberOfArrays()):
                a = info.GetArray(i)
                if a.GetNumberOfComponents() >= 3 and a.GetName() not in point_names:
                    gvec = a.GetName()
                    break
        except Exception:
            pass

    if gvec is None:
        print("  WARNING: no gradient array found; arrays=%s" % gnames, file=sys.stderr)
    else:
        calc = Calculator(Input=grad)
        calc.ResultArrayName = "schlieren"
        calc.Function = "mag(%s)" % gvec
        UpdatePipeline(proxy=calc)
        _hide_all(view)
        d2 = Show(calc, view)
        d2.Representation = "Surface"
        ColorBy(d2, ("POINTS", "schlieren"))
        d2.RescaleTransferFunctionToDataRange(True, False)
        d2.SetScalarBarVisibility(view, True)
        _apply_lut("schlieren", ["Grayscale", "X Ray", "Black-Body Radiation"], invert=True)
        _top_down_camera(view, bounds)
        _save(
            view,
            os.path.join(outdir, "%s_schlieren.png" % prefix),
            args.width,
            args.height,
            args.res,
        )
        try:
            Delete(calc)
            Delete(grad)
        except Exception:
            pass

    # ---- Sensor / AV ----
    print("Sensor / AV field…")
    field = None
    assoc = "POINTS"
    if "sensor" in point_names:
        field, assoc = "sensor", "POINTS"
    elif "sensor" in cell_names:
        field, assoc = "sensor", "CELLS"
    elif "av" in point_names:
        field, assoc = "av", "POINTS"
    elif "av" in cell_names:
        field, assoc = "av", "CELLS"

    if field is None:
        print(
            "  WARNING: no sensor/av arrays in VTU. Re-run with diagnostics:\n"
            "    write_vtu_high_order_with_capturing(...) or --vtk-diagnostics"
        )
    else:
        _hide_all(view)
        d3 = Show(reader, view)
        d3.Representation = "Surface"
        ColorBy(d3, (assoc, field))
        d3.RescaleTransferFunctionToDataRange(True, False)
        d3.SetScalarBarVisibility(view, True)
        _apply_lut(field, ["Plasma (matplotlib)", "Cool to Warm", "Rainbow Desaturated"])
        _top_down_camera(view, bounds)
        _save(
            view,
            os.path.join(outdir, "%s_sensor.png" % prefix),
            args.width,
            args.height,
            args.res,
        )

    # ---- Line-outs ----
    print("Line-outs…")
    lines = _lineout_presets(case, bounds)
    _export_lineouts(
        reader,
        lines,
        os.path.join(outdir, "%s_lineout.csv" % prefix),
        os.path.join(outdir, "%s_lineout.png" % prefix),
        args.width,
        args.height,
        args.res,
    )

    print("Done. Outputs under", outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
