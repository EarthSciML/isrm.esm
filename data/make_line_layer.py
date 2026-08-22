#!/usr/bin/env python3
"""Build the example LINE emission layer `isrm_line.esm` reads.

The `.esm` document states physics, not data preparation: it reads a shapefile
through EarthSciIO's `shapefile` reader and allocates each road's emissions to
the ISRM source cells by the fraction of the road's length that falls in each
cell.  THIS script is the other half - the step a user does in Python before
pointing a document at the result - and it is deliberately ordinary: fetch a
public road file, keep the roads you care about, thin them to the grid's own
resolution, cut them into segments, attach an emission column, write a
shapefile.

    python data/make_line_layer.py                    # Illinois interstates
    python data/make_line_layer.py --state 06         # California
    python data/make_line_layer.py --simplify-m 500   # coarser shape
    python data/make_line_layer.py --rate 2.5 --out /tmp/mylayer.zip

What it writes
--------------
One `.zip` holding a four-file shapefile set - `.shp` (geometry), `.shx`
(index), `.dbf` (attributes), `.prj` (the CRS as WKT), exactly as
`make_polygon_layer.py` writes its counterpart.

ONE ROW PER SEGMENT
-------------------
Every row is a two-vertex polyline, and the rows of one road carry that road's
`LINEID` and its whole `EMIS`.  This is the layer's one non-obvious choice, and
it is forced by the engine rather than chosen for convenience.

A road is a polyline with many vertices, so the natural document would range an
emission-binning aggregate over (cell, road, segment) and contract the last two
away.  That is legal ESM - the spec's aggregate contracts a SET of index
symbols - but the projection-pushdown rewrite that makes this model tractable
refuses it: `_pd_detect_binning` requires the binning aggregate to declare
EXACTLY TWO ranges, one of which is the output cell set.  A three-range
aggregate is not recognised, no support set is derived, and the 33 GB
source-receptor slab is fetched whole.  So the segment, not the road, has to be
the record - and the exploding has to happen here, because the reader delivers
what the file holds.

What the document still does for itself: every length in the model - each
segment's, each road's total, and the clipped length inside a cell - is
computed by `isrm_line.esm` from the projected geometry.  This script hands it
coordinates and a parent id, not lengths.  `LEN_KM` below is written for a
reader of the layer and is deliberately not what the document divides by; see
the `line_len` note in the document.

Columns
-------
`LINEID`   parent road id: rows sharing it are one road, and the document sums
           their lengths to get the denominator of the length share.
`SEGID`    the segment's position along its road, from 1.  Not read by the
           document; it makes a row identifiable when looking at the layer.
`FULLNAME` the road's name from TIGER (`I- 55`), replicated across its segments.
`EMIS`     the PARENT ROAD's total emission, short tons per year, replicated
           across its segments.  Replicated and not divided up, because
           dividing it here would be the allocation the document exists to do.
`LEN_KM`   the parent road's geodesic length in km, replicated.  Written for a
           human reading the layer; the document computes its own planar length
           instead, and the ~0.2% the two differ by is the map projection.

The emission is DECLARED, not measured
--------------------------------------
`EMIS = rate * LEN_KM` short tons per year: a uniform per-kilometre line
source, at a rate this script's `--rate` sets.  It is an EXAMPLE quantity and
the document says so - the point of the line path is the geometry (how a line's
mass reaches the cells it crosses), not the inventory.  `isrm_point.esm` is
where a real, sourced inventory lives.

The shape is THINNED, and that is declared too
----------------------------------------------
TIGER road centrelines carry a vertex every ~50 m.  The ISRM grid's finest cell
is 1 km, so that detail cannot reach the model - it only makes the record axis
50x longer than the answer needs.  Each road is therefore simplified by
Douglas-Peucker at `--simplify-m` (default 200 m, a fifth of the finest cell),
which is a bounded, stated approximation: no point of the kept polyline is more
than the tolerance from the original.  `EMIS` is derived from the length of the
SIMPLIFIED road, so the layer is self-consistent - the emission and the geometry
describe the same object - and the log prints what the thinning cost.
"""

from __future__ import annotations

import argparse
import io
import math
import os
import pathlib
import sys
import urllib.request
import zipfile

# The Census TIGER/Line primary and secondary roads for one state - the
# published road centrelines, already a zipped shapefile set.  Pinned to the
# 2020 vintage so a rebuild is reproducible.
SOURCE_TEMPLATE = (
    "https://www2.census.gov/geo/tiger/TIGER2020/PRISECROADS/"
    "tl_2020_{state}_prisecroads.zip"
)
DEFAULT_STATE = "17"  # Illinois
DEFAULT_RATE = 1.0  # short tons per year per kilometre of road
DEFAULT_SIMPLIFY = 200.0  # metres; a fifth of the ISRM grid's finest cell
DEFAULT_MTFCC = "S1100"  # primary road (interstate / US highway)
DEFAULT_RTTYP = "I"  # route type: Interstate

R_EARTH = 6371008.8  # m, the IUGG mean radius - for lengths only, never for the model


def fetch(url: str, cache: pathlib.Path) -> bytes:
    """The source zip's bytes, downloaded once into `cache`."""
    if cache.is_file():
        return cache.read_bytes()
    print(f"  downloading {url}")
    with urllib.request.urlopen(url, timeout=600) as resp:
        data = resp.read()
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_bytes(data)
    return data


def local_metres(pts: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """`pts` (lon, lat degrees) in a local equirectangular frame, in metres.

    Used only to give `--simplify-m` and the length columns a unit.  It is not
    the model's projection: the document reprojects the lon/lat this script
    writes with its own Lambert conformal conic, and never sees this frame.
    """
    lat0 = math.radians(sum(p[1] for p in pts) / len(pts))
    k = math.radians(1.0) * R_EARTH
    return [(x * k * math.cos(lat0), y * k) for x, y in pts]


def geodesic_length(pts: list[tuple[float, float]]) -> float:
    """Length of the lon/lat polyline `pts`, in metres."""
    total = 0.0
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        lat = math.radians((y0 + y1) / 2.0)
        dx = math.radians(x1 - x0) * math.cos(lat) * R_EARTH
        dy = math.radians(y1 - y0) * R_EARTH
        total += math.hypot(dx, dy)
    return total


def simplify(pts: list[tuple[float, float]], tol_m: float) -> list[tuple[float, float]]:
    """Douglas-Peucker over `pts` (lon, lat) at a tolerance in METRES.

    Endpoints are always kept, and no dropped vertex is further than `tol_m`
    from the surviving polyline.  Iterative rather than recursive: a TIGER road
    can carry 1,700 vertices and Python's recursion limit is not the place to
    discover that.
    """
    if tol_m <= 0 or len(pts) < 3:
        return list(pts)
    m = local_metres(pts)
    keep = [False] * len(pts)
    keep[0] = keep[-1] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        a, b = stack.pop()
        if b <= a + 1:
            continue
        x0, y0 = m[a]
        x1, y1 = m[b]
        dx, dy = x1 - x0, y1 - y0
        norm = math.hypot(dx, dy)
        worst, worst_i = -1.0, -1
        for i in range(a + 1, b):
            px, py = m[i]
            if norm > 0:  # distance from the chord, or from the point it collapsed to
                d = abs(dx * (y0 - py) - dy * (x0 - px)) / norm
            else:
                d = math.hypot(px - x0, py - y0)
            if d > worst:
                worst, worst_i = d, i
        if worst > tol_m:
            keep[worst_i] = True
            stack.append((a, worst_i))
            stack.append((worst_i, b))
    return [p for p, k in zip(pts, keep) if k]


def read_parts(source: bytes, mtfcc: str, rttyp: str):
    """The source layer's matching PARTS: `(name, [(lon, lat), ...])` each."""
    import shapefile  # pyshp - the same decoder EarthSciIO's reader uses

    with zipfile.ZipFile(io.BytesIO(source)) as zf:
        names = [n for n in zf.namelist() if not n.endswith("/")]
        shp = next(n for n in names if n.lower().endswith(".shp"))
        stem = shp[: -len(".shp")]
        blobs = {
            ext: zf.read(f"{stem}.{ext}")
            for ext in ("shp", "shx", "dbf", "prj")
            if f"{stem}.{ext}" in names
        }

    reader = shapefile.Reader(
        shp=io.BytesIO(blobs["shp"]),
        shx=io.BytesIO(blobs["shx"]),
        dbf=io.BytesIO(blobs["dbf"]),
    )
    parts = []
    for sr in reader.iterShapeRecords():
        rec = sr.record.as_dict()
        if mtfcc and str(rec.get("MTFCC", "")) != mtfcc:
            continue
        if rttyp and str(rec.get("RTTYP", "")) != rttyp:
            continue
        pts = [(float(x), float(y)) for x, y in sr.shape.points]
        bounds = list(sr.shape.parts) or [0]
        name = str(rec.get("FULLNAME", "") or "").strip()
        for i, start in enumerate(bounds):
            stop = bounds[i + 1] if i + 1 < len(bounds) else len(pts)
            if stop - start >= 2:
                parts.append((name, pts[start:stop]))
    reader.close()
    return parts, blobs.get("prj")


def build(parts, prj, rate: float, tol_m: float):
    """The output shapefile's members, plus a summary for the log."""
    import shapefile

    shp_io, shx_io, dbf_io = io.BytesIO(), io.BytesIO(), io.BytesIO()
    raw_len = 0.0
    kept_len = 0.0
    n_seg = 0
    raw_vertices = 0
    kept_vertices = 0
    with shapefile.Writer(shp=shp_io, shx=shx_io, dbf=dbf_io, shapeType=shapefile.POLYLINE) as w:
        w.field("LINEID", "N", 9, 0)
        w.field("SEGID", "N", 9, 0)
        w.field("FULLNAME", "C", 40)
        w.field("EMIS", "N", 18, 6)
        w.field("LEN_KM", "N", 18, 6)
        line_id = 0
        for name, pts in parts:
            raw_len += geodesic_length(pts)
            raw_vertices += len(pts)
            thin = simplify(pts, tol_m)
            # A zero-length segment carries no mass and would divide by zero in
            # any per-segment normalisation, so it never reaches the layer.
            segs = [(a, b) for a, b in zip(thin, thin[1:]) if a != b]
            if not segs:
                continue
            length = geodesic_length(thin)
            if length <= 0.0:
                continue
            line_id += 1
            kept_len += length
            kept_vertices += len(thin)
            emis = rate * length / 1000.0
            for k, (a, b) in enumerate(segs, start=1):
                w.line([[list(a), list(b)]])
                w.record(line_id, k, name, emis, length / 1000.0)
                n_seg += 1

    if n_seg == 0:
        raise SystemExit("no road survived the filters - nothing to write")

    out = {
        "emis_lines.shp": shp_io.getvalue(),
        "emis_lines.shx": shx_io.getvalue(),
        "emis_lines.dbf": dbf_io.getvalue(),
    }
    if prj is not None:  # carry the source CRS forward verbatim
        out["emis_lines.prj"] = prj
    summary = dict(lines=line_id, segments=n_seg, raw_len=raw_len, kept_len=kept_len,
                   raw_vertices=raw_vertices, kept_vertices=kept_vertices,
                   emis=rate * kept_len / 1000.0)
    return out, summary


def write_zip(path: pathlib.Path, members: dict[str, bytes]) -> None:
    """Write the shapefile set as one zip, byte-deterministically.

    Fixed timestamps + STORED entries mean rebuilding the same layer produces
    the same bytes, so the content-addressed cache sees one blob rather than a
    new one per rebuild.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w") as zf:
        for name in sorted(members):
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.create_system = 3
            info.external_attr = 0o644 << 16
            info.compress_type = zipfile.ZIP_STORED
            zf.writestr(info, members[name])


def main(argv: list[str]) -> int:
    here = pathlib.Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--state", default=DEFAULT_STATE,
                    help="two-digit state FIPS to build (default 17 = Illinois)")
    ap.add_argument("--rate", type=float, default=DEFAULT_RATE,
                    help="short tons per year per KM of road (default 1.0)")
    ap.add_argument("--simplify-m", type=float, default=DEFAULT_SIMPLIFY,
                    help="Douglas-Peucker tolerance in metres, 0 to keep every "
                         "vertex (default 200)")
    ap.add_argument("--mtfcc", default=DEFAULT_MTFCC,
                    help="TIGER feature class to keep, '' for all "
                         "(default S1100 = primary road)")
    ap.add_argument("--rttyp", default=DEFAULT_RTTYP,
                    help="TIGER route type to keep, '' for all (default I = interstate)")
    ap.add_argument("--source", default=None,
                    help="road zip to read (URL or local path; default is the "
                         "TIGER2020 PRISECROADS file for --state)")
    ap.add_argument("--out", default=None,
                    help="output zip (default data/line_emissions_<state>.zip)")
    args = ap.parse_args(argv[1:])

    try:
        import shapefile  # noqa: F401
    except ImportError:
        print("this script needs pyshp:  pip install pyshp", file=sys.stderr)
        return 2

    src = args.source or SOURCE_TEMPLATE.format(state=args.state)
    out = pathlib.Path(args.out) if args.out else (
        here / f"line_emissions_{args.state}.zip")
    print(f"building {out.name}")
    if os.path.isfile(src):
        source = pathlib.Path(src).read_bytes()
    else:
        source = fetch(src, here / f"tl_2020_{args.state}_prisecroads.zip")

    parts, prj = read_parts(source, args.mtfcc, args.rttyp)
    if not parts:
        raise SystemExit(f"no part matches MTFCC {args.mtfcc!r} / RTTYP {args.rttyp!r}")
    members, s = build(parts, prj, args.rate, args.simplify_m)
    write_zip(out, members)

    lost = 100.0 * (s["raw_len"] - s["kept_len"]) / s["raw_len"]
    print(f"  {s['lines']} roads -> {s['segments']} segments "
          f"(one row each, two vertices each)")
    print(f"  Douglas-Peucker at {args.simplify_m:g} m: "
          f"{s['raw_vertices']} vertices -> {s['kept_vertices']}, "
          f"length {s['raw_len'] / 1000:,.1f} -> {s['kept_len'] / 1000:,.1f} km "
          f"({lost:.3f}% shorter)")
    print(f"  {s['kept_len'] / 1000:,.1f} km of road -> {s['emis']:,.1f} short tons/yr "
          f"at {args.rate} t/yr/km")
    print(f"  wrote {out}  ({out.stat().st_size:,} bytes)")
    print("  columns: LINEID, SEGID, FULLNAME, EMIS (parent road's short tons/yr), "
          "LEN_KM (parent road's km)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
