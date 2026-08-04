#!/usr/bin/env python3
"""
Publication-style 2D figures from FRForge high-order VTU (ParaView / pvpython).

Produces (for a given --prefix):
  {prefix}_schlieren.png   — inverted grayscale |∇ρ|
  {prefix}_pressure.png    — pressure field
  {prefix}_sensor.png      — modal sensor σ (or av if sensor missing)
  {prefix}_lineout.png     — density (and pressure) line cuts
  {prefix}_lineout.csv     — numeric samples for the cuts

Requirements:
  - ParaView with pvpython (or `paraview` Python package matching your install)
  - VTU from FRForge write_vtu_high_order[_with_capturing]

Examples:
  pvpython scripts/docs/paraview/plot_2d_publication.py \\
    --vtu results/docs_vtu/riemann_cfg6_persson_av.vtu \\
    --outdir results/docs_figures \\
    --prefix riemann_cfg6_baseline \\
    --case riemann

  pvpython scripts/docs/paraview/plot_2d_publication.py \\
    --vtu results/docs_vtu/double_mach_persson_av.vtu \\
    --outdir results/docs_figures \\
    --prefix double_mach_baseline \\
    --case dmr
"""

from __future__ import print_function

import argparse
import csv
import os
import sys


def _import_paraview():
    try:
        from paraview.simple import (
            XMLUnstructuredGridReader,
            Gradient,
            Calculator,
            Contour,
            PlotOverLine,
            GetActiveViewOrCreate,
            GetColorTransferFunction,
            GetOpacityTransferFunction,
            GetAnimationScene,
            ColorBy,
            Hide,
            Show,
            SaveScreenshot,
            UpdatePipeline,
            Delete,
            Render,
            ResetCamera,
            CreateLayout,
        )
        import paraview.servermanager as sm
        return True
    except Exception as e:
        print(
            "ERROR: Could not import ParaView Python bindings.\n"
            "  Run this script with `pvpython` from your ParaView install, e.g.:\n"
            "    /Applications/ParaView-*.app/Contents/bin/pvpython ...\n"
            "    or:  pvpython scripts/docs/paraview/plot_2d_publication.py ...\n"
            "Import error: %s" % e,
            file=sys.stderr,
        )
        return False


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


def _point_array_names(reader):
    """Return list of point data array names."""
    try:
        info = reader.GetPointDataInformation()
        names = []
        for i in range(info.GetNumberOfArrays()):
            names.append(info.GetArray(i).GetName())
        return names
    except Exception:
        return []


def _cell_array_names(reader):
    try:
        info = reader.GetCellDataInformation()
        names = []
        for i in range(info.GetNumberOfArrays()):
            names.append(info.GetArray(i).GetName())
        return names
    except Exception:
        return []


def _setup_view(size=(1400, 1200)):
    from paraview.simple import GetActiveViewOrCreate

    view = GetActiveViewOrCreate("RenderView")
    view.ViewSize = list(size)
    view.Background = [1.0, 1.0, 1.0]
    view.OrientationAxesVisibility = 0
    view.CameraParallelProjection = 1
    return view


def _hide_all(view):
    from paraview.simple import GetSources, Hide

    for _, src in GetSources().items():
        try:
            Hide(src, view)
        except Exception:
            pass


def _color_by_point(display, view, name, cmap="Cool to Warm", invert=False, log=False):
    from paraview.simple import (
        ColorBy,
        GetColorTransferFunction,
        GetOpacityTransferFunction,
    )

    ColorBy(display, ("POINTS", name))
    display.RescaleTransferFunctionToDataRange(True, False)
    display.SetScalarBarVisibility(view, True)
    lut = GetColorTransferFunction(name)
    try:
        lut.ApplyPreset(cmap, True)
    except Exception:
        pass
    if invert:
        try:
            lut.InvertTransferFunction()
        except Exception:
            # manual reverse: swap RGB points if needed — Invert is available in modern PV
            pass
    if log:
        try:
            lut.MapControlPointsToLogSpace()
            lut.UseLogScale = 1
        except Exception:
            pass
    # clean scalar bar
    try:
        sb = display.LookupTable
        # ScalarBar is accessed via view representations
        for proxy in view.Representations:
            pass
    except Exception:
        pass
    return lut


def _save(view, path, width, height, res=2):
    from paraview.simple import SaveScreenshot, Render

    Render(view)
    SaveScreenshot(
        path,
        view,
        ImageResolution=[width * res, height * res],
        OverrideColorPalette="WhiteBackground",
        TransparentBackground=0,
        CompressionLevel="5",
    )
    print("  wrote", path)


def _bounds(reader):
    from paraview.simple import UpdatePipeline

    UpdatePipeline()
    d = reader.GetDataInformation().GetBounds()
    # (xmin,xmax,ymin,ymax,zmin,zmax)
    return d


def _lineout_presets(case, bounds):
    xmin, xmax, ymin, ymax = bounds[0], bounds[1], bounds[2], bounds[3]
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)
    if case == "dmr":
        # horizontal through lower flow; vertical mid-domain
        lines = [
            ("y_mid", (xmin, cy, 0.0), (xmax, cy, 0.0)),
            ("x_shock", (cx, ymin, 0.0), (cx, ymax, 0.0)),
        ]
    else:
        # Riemann unit square defaults
        lines = [
            ("y_0.5", (xmin, cy, 0.0), (xmax, cy, 0.0)),
            ("x_0.5", (cx, ymin, 0.0), (cx, ymax, 0.0)),
        ]
    return lines


def _export_lineouts(reader, lines, csv_path, png_path, view, width, height, res):
    from paraview.simple import (
        PlotOverLine,
        Show,
        Hide,
        UpdatePipeline,
        CreateView,
        GetActiveViewOrCreate,
        SaveScreenshot,
        Render,
        Delete,
    )

    rows = []
    headers = ["line", "s", "x", "y", "rho", "p"]
    # Use first line for PNG XY chart if available
    first_plot = None
    for name, p0, p1 in lines:
        pol = PlotOverLine(Input=reader)
        pol.Point1 = list(p0)
        pol.Point2 = list(p1)
        # denser sampling
        try:
            pol.Resolution = 500
        except Exception:
            pass
        UpdatePipeline()
        # fetch data via servermanager
        try:
            data = pol.GetClientSideObject().GetOutput()
            # vtkTable or vtkPolyData
            import vtk

            pdata = pol.GetClientSideObject().GetOutputDataObject(0)
            # Prefer Fetch
            from paraview.simple import servermanager as sm

            fetcher = sm.Fetch(pol)
            # vtkPolyData with point data
            n = fetcher.GetNumberOfPoints()
            pd = fetcher.GetPointData()
            arr_rho = pd.GetArray("rho")
            arr_p = pd.GetArray("p")
            for i in range(n):
                pt = fetcher.GetPoint(i)
                # arc length param ~ i/(n-1)
                s = float(i) / max(n - 1, 1)
                rho = arr_rho.GetValue(i) if arr_rho is not None else float("nan")
                pr = arr_p.GetValue(i) if arr_p is not None else float("nan")
                rows.append([name, s, pt[0], pt[1], rho, pr])
        except Exception as e:
            print("  WARNING: line-out fetch failed for", name, ":", e)

        if first_plot is None:
            first_plot = pol
        else:
            try:
                Delete(pol)
            except Exception:
                pass

    # CSV
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(headers)
        w.writerows(rows)
    print("  wrote", csv_path)

    # Simple line chart view for first cut if possible
    if first_plot is not None:
        try:
            chart = CreateView("XYChartView")
            chart.ViewSize = [width, height]
            disp = Show(first_plot, chart)
            try:
                disp.SeriesVisibility = ["rho", "p"]
            except Exception:
                pass
            Render(chart)
            SaveScreenshot(
                png_path,
                chart,
                ImageResolution=[width * res, height * res],
                OverrideColorPalette="WhiteBackground",
            )
            print("  wrote", png_path)
            Delete(chart)
        except Exception as e:
            print("  WARNING: could not save lineout PNG (CSV is still written):", e)
            # fallback empty note
            if not os.path.isfile(png_path):
                print("  (skip PNG)")


def main(argv=None):
    if not _import_paraview():
        return 2

    from paraview.simple import (
        XMLUnstructuredGridReader,
        GradientOfUnstructuredDataSet,
        Calculator,
        Show,
        Hide,
        UpdatePipeline,
        ColorBy,
        GetColorTransferFunction,
        Render,
        ResetCamera,
        Delete,
        GetActiveViewOrCreate,
    )

    # Gradient filter name differs slightly across versions
    try:
        from paraview.simple import Gradient
        gradient_filter = Gradient
    except Exception:
        gradient_filter = GradientOfUnstructuredDataSet

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
    UpdatePipeline()
    point_names = _point_array_names(reader)
    cell_names = _cell_array_names(reader)
    print("  PointData:", point_names)
    print("  CellData:", cell_names)

    if "rho" not in point_names:
        print("ERROR: PointData 'rho' required for Schlieren/pressure plots.", file=sys.stderr)
        return 1

    # case preset
    case = args.case
    if case == "auto":
        case = "dmr" if "double_mach" in os.path.basename(vtu).lower() or "dmr" in prefix.lower() else "riemann"

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
        lut = GetColorTransferFunction("p")
        try:
            lut.ApplyPreset("Cool to Warm", True)
        except Exception:
            try:
                lut.ApplyPreset("Viridis (matplotlib)", True)
            except Exception:
                pass
    ResetCamera(view)
    # slightly pad camera
    try:
        view.CameraParallelScale = view.CameraParallelScale * 1.05
    except Exception:
        pass
    _save(
        view,
        os.path.join(outdir, "%s_pressure.png" % prefix),
        args.width,
        args.height,
        args.res,
    )

    # ---- Schlieren |∇ρ| inverted grayscale ----
    print("Numerical Schlieren…")
    _hide_all(view)
    grad = gradient_filter(Input=reader)
    # attribute selection differs by filter version
    try:
        grad.ScalarArray = ["POINTS", "rho"]
    except Exception:
        try:
            grad.SelectInputScalars = ["POINTS", "rho"]
        except Exception:
            pass
    try:
        grad.ComputeGradient = 1
    except Exception:
        pass
    UpdatePipeline()

    # Magnitude via Calculator if Result array is vector
    calc = Calculator(Input=grad)
    calc.ResultArrayName = "schlieren"
    # Try common gradient array names
    gnames = _point_array_names(grad)
    gvec = None
    for cand in ("Gradients", "Gradient", "rhoGradient", "Result"):
        if cand in gnames:
            gvec = cand
            break
    if gvec is None and gnames:
        gvec = gnames[-1]
    if gvec is None:
        print("  WARNING: no gradient array; skipping Schlieren", file=sys.stderr)
    else:
        calc.Function = "mag(%s)" % gvec
        UpdatePipeline()
        _hide_all(view)
        d2 = Show(calc, view)
        d2.Representation = "Surface"
        ColorBy(d2, ("POINTS", "schlieren"))
        d2.RescaleTransferFunctionToDataRange(True, False)
        d2.SetScalarBarVisibility(view, True)
        lut = GetColorTransferFunction("schlieren")
        try:
            lut.ApplyPreset("Grayscale", True)
        except Exception:
            try:
                lut.ApplyPreset("X Ray", True)
            except Exception:
                pass
        # Invert so high gradient = black
        try:
            lut.InvertTransferFunction()
        except Exception:
            pass
        # Prefer linear scale; clamp high outliers via custom range if needed
        ResetCamera(view)
        try:
            view.CameraParallelScale = view.CameraParallelScale * 1.05
        except Exception:
            pass
        _save(
            view,
            os.path.join(outdir, "%s_schlieren.png" % prefix),
            args.width,
            args.height,
            args.res,
        )

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
        lut = GetColorTransferFunction(field)
        try:
            lut.ApplyPreset("Plasma (matplotlib)", True)
        except Exception:
            try:
                lut.ApplyPreset("Cool to Warm", True)
            except Exception:
                pass
        ResetCamera(view)
        try:
            view.CameraParallelScale = view.CameraParallelScale * 1.05
        except Exception:
            pass
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
        view,
        args.width,
        args.height,
        args.res,
    )

    print("Done. Outputs under", outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
