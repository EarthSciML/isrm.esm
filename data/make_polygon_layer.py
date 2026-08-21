#!/usr/bin/env python3
"""Build the example POLYGON emission layer `isrm_polygon.esm` reads.

The `.esm` document states physics, not data preparation: it reads a shapefile
through EarthSciIO's `shapefile` reader and allocates each polygon's emissions
to the ISRM source cells by area overlap. THIS script is the other half — the
step a user does in Python before pointing a document at the result — and it is
deliberately ordinary: fetch a public boundary file, keep the polygons you care
about, attach an emission column, write a shapefile.

    python data/make_polygon_layer.py                 # Illinois, 1 t/yr/km^2
    python data/make_polygon_layer.py --state 06      # California (islands!)
    python data/make_polygon_layer.py --rate 2.5 --out /tmp/mylayer.zip

What it writes
--------------
One `.zip` holding a four-file shapefile set — `.shp` (geometry), `.shx`
(index), `.dbf` (attributes), `.prj` (the CRS as WKT). The zip is the fetchable
form because the content-addressed cache holds ONE blob per URL, and it is what
the reader's `member` option names a layer inside.

Columns: `GEOID` (5-digit county FIPS, kept as text — the reader's
`numeric_columns` option turns it into a number when a model wants one), `NAME`,
`ALAND` (land area, m^2, straight from the source file) and `EMIS`, the emission
this script derives. Geometry stays in the source file's own CRS (NAD83
geographic, EPSG:4269 — the same lon/lat datum the FF10 point inventory uses),
so a document reprojects it with the Lambert-conformal template it already owns.

The emission is DECLARED, not measured
--------------------------------------
`EMIS = rate * ALAND / 1e6` short tons per year: a uniform per-square-kilometre
source, at a rate this script's `--rate` sets. It is an EXAMPLE quantity and the
document says so — the point of the polygon path is the geometry (how a polygon's
mass reaches the grid cells it covers), not the inventory. `isrm_point.esm` is
where a real, sourced inventory lives.
"""

from __future__ import annotations

import argparse
import io
import os
import pathlib
import sys
import urllib.request
import zipfile

# The US Census cartographic boundary counties, 1:20,000,000 — the smallest
# national county layer the Census publishes (~890 KB), already a zipped
# shapefile set. Pinned to the 2020 vintage so a rebuild is reproducible.
SOURCE_URL = (
    "https://www2.census.gov/geo/tiger/GENZ2020/shp/cb_2020_us_county_20m.zip"
)
DEFAULT_STATE = "17"  # Illinois
DEFAULT_RATE = 1.0  # short tons per year per square kilometre of land


def fetch(url: str, cache: pathlib.Path) -> bytes:
    """The source zip's bytes, downloaded once into `cache`."""
    if cache.is_file():
        return cache.read_bytes()
    print(f"  downloading {url}")
    with urllib.request.urlopen(url, timeout=300) as resp:
        data = resp.read()
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_bytes(data)
    return data


def build(source: bytes, state: str, rate: float) -> tuple[dict[str, bytes], list]:
    """The output shapefile's members, plus a per-record summary for the log."""
    import shapefile  # pyshp — the same decoder EarthSciIO's reader uses

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
    shp_io, shx_io, dbf_io = io.BytesIO(), io.BytesIO(), io.BytesIO()
    summary = []
    with shapefile.Writer(shp=shp_io, shx=shx_io, dbf=dbf_io) as w:
        w.field("GEOID", "C", 5)
        w.field("NAME", "C", 32)
        w.field("ALAND", "N", 18, 0)
        w.field("EMIS", "N", 18, 6)
        for sr in reader.iterShapeRecords():
            rec = sr.record.as_dict()
            if str(rec["STATEFP"]) != state:
                continue
            aland = float(rec["ALAND"])
            emis = rate * aland / 1e6
            w.shape(sr.shape)
            w.record(str(rec["GEOID"]), str(rec["NAME"]), aland, emis)
            summary.append((rec["GEOID"], rec["NAME"], len(sr.shape.parts) or 1,
                            len(sr.shape.points), aland, emis))
    reader.close()
    if not summary:
        raise SystemExit(f"no county has STATEFP {state!r} — nothing to write")

    out = {
        "emis_polygons.shp": shp_io.getvalue(),
        "emis_polygons.shx": shx_io.getvalue(),
        "emis_polygons.dbf": dbf_io.getvalue(),
    }
    if "prj" in blobs:  # carry the source CRS forward verbatim
        out["emis_polygons.prj"] = blobs["prj"]
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
                    help="two-digit state FIPS to keep (default 17 = Illinois)")
    ap.add_argument("--rate", type=float, default=DEFAULT_RATE,
                    help="short tons per year per km^2 of LAND area (default 1.0)")
    ap.add_argument("--source", default=SOURCE_URL,
                    help="county boundary zip to read (URL or local path)")
    ap.add_argument("--out", default=None,
                    help="output zip (default data/polygon_emissions_<state>.zip)")
    args = ap.parse_args(argv[1:])

    try:
        import shapefile  # noqa: F401
    except ImportError:
        print("this script needs pyshp:  pip install pyshp", file=sys.stderr)
        return 2

    out = pathlib.Path(args.out) if args.out else (
        here / f"polygon_emissions_{args.state}.zip")
    print(f"building {out.name}")
    if os.path.isfile(args.source):
        source = pathlib.Path(args.source).read_bytes()
    else:
        source = fetch(args.source, here / "cb_2020_us_county_20m.zip")

    members, summary = build(source, args.state, args.rate)
    write_zip(out, members)

    parts = sum(s[2] for s in summary)
    print(f"  {len(summary)} polygons ({parts} ring(s)), "
          f"{max(s[3] for s in summary)} vertices in the largest record")
    print(f"  land area {sum(s[4] for s in summary) / 1e6:,.0f} km^2 "
          f"-> {sum(s[5] for s in summary):,.1f} short tons/yr at "
          f"{args.rate} t/yr/km^2")
    print(f"  wrote {out}  ({out.stat().st_size:,} bytes)")
    print("  columns: GEOID (text FIPS), NAME, ALAND (m^2), EMIS (short tons/yr)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
