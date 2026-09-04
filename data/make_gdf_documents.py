#!/usr/bin/env python3
"""Generate the GeoDataFrame-driven ISRM documents from the published ones.

    python data/make_gdf_documents.py            # write them
    python data/make_gdf_documents.py --check    # fail if they are out of date

`isrm_point.esm`, `isrm_polygon.esm` and `isrm_line.esm` are the cross-language
contract: three shims drive them in-process, and their totals are checked
against published InMAP results. Nothing here edits them. What this script
writes is a SIBLING SET for a different caller — `run-api/isrm_esm.ipynb`, whose
one function takes a `geopandas` GeoDataFrame, uploads it as a shapefile and
runs it through EarthSciLab's API:

    isrm_gdf_point.esm      point sources, ASME plume rise across three layers
    isrm_gdf_point_flat.esm the same without stack parameters: layer 0 only
    isrm_gdf_polygon.esm    area sources, allocated by area overlap
    isrm_gdf_line.esm       line sources, allocated by within-cell length

Three things differ from their parents, and each one is forced by the caller
rather than chosen for tidiness.

ONE SHAPEFILE, FIVE POLLUTANT COLUMNS
-------------------------------------
A GeoDataFrame is WIDE: one row per source, one column per pollutant. FF10 is
LONG — one row per (stack, pollutant) — which is why `isrm_point.esm` carries a
77-entry pollutant-code map and five `is_VOC`/`is_NOx`/... membership masks to
recover the pathway a row feeds. Given a column per pollutant those masks say
nothing: `E_PM25_L0`'s `emis_annual[r] * is_PM25[r]` becomes `emis_PM25[r]`, and
the code map, the masks and the `pollutant` variable all leave the document.
This set is SMALLER than its parent, not larger.

The area documents go the other way. Both parents wire exactly one pathway, and
`isrm_polygon.esm`'s own `TotalPM25` says a five-pathway version "adds terms
here and nowhere else" — so the five-pathway layer is lifted off
`isrm_point.esm` (5 `SR_*_L0`, 5 `conc_*`, 5 `pm_*`, a five-term `TotalPM25`)
and set on top of each parent's own binning equation, which is the only thing
that distinguishes area from line from point.

SI UNITS
--------
The caller hands over kg/yr, metres, kelvin and m/s, so these documents declare
those. `fact` — the emission-to-µg/s conversion applied once to the summed
concentration — drops from 28766.639396245562 (short_ton/yr) to
1e9/3.1536e7 = 31.709791983764585 (kg/yr): exactly the parent's constant with
the 907.18474 kg/short-ton factor divided out, so nothing else moves. The point
document already computed `stk_height` (m) from `stkhgt` (ft); now it reads the
metres directly and those four conversion equations go.

THE ROAD -> SEGMENT SPLIT MOVES TO THE CALLER
---------------------------------------------
`isrm_line.esm` receives whole roads and apportions each road's emission across
its segments inside the document: `line_key` groups the segments, `road_len`
sums each road's length with an aggregate ranging over `emis_records` TWICE,
and `seg_emis` divides. That O(N_REC^2) group-total is only there because the
reader delivers what the file holds. A caller holding a GeoDataFrame explodes
the road and apportions by length before it ever writes the file, so these
documents read per-segment emissions directly and `line_key`, `road_len` and
`seg_emis` are gone with the quadratic aggregate.
"""

from __future__ import annotations

import argparse
import copy
import json
import pathlib
import sys

import earthsci_ast

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent

# Each emission column and the SR pathway it feeds. The column names are the
# caller's API; the pathway names are the ISRM store's own array names.
PATHWAYS = {
    "PM25": "PrimaryPM25",
    "VOC": "SOA",
    "NOx": "pNO3",
    "NH3": "pNH4",
    "SOx": "pSO4",
}

# kg/yr -> ug/s. The parent's 28766.639396245562 is this times 907.18474 kg per
# short ton, so a document that reads kg/yr is numerically the parent.
FACT_SI = 1e9 / 3.1536e7

# The stack parameters ASME plume rise needs, in the units the plume-rise
# expressions already work in, mapped to the shapefile column that carries each.
STACK_COLUMNS = {
    "stk_height": ("STKHGT", "m", "Stack height above ground."),
    "stk_diam": ("STKDIAM", "m", "Stack inner diameter at the exit."),
    "stk_temp": ("STKTEMP", "K", "Exit gas temperature."),
    "stk_vel": ("STKVEL", "m/s", "Exit gas velocity."),
}

SOURCE = "Emis"          # every generated document names its emission loader this
MEMBER = "emis.shp"      # and expects this .shp inside the uploaded zip

READER_NOTE = (
    "The reader is READER-ONLY: it decodes rings and .dbf columns and declares the "
    "file's CRS, and it reprojects nothing. `member` names the .shp inside the zip "
    "(the .dbf/.shx/.prj sidecars are the same stem) and never enters the cache key. "
    "`nvert_max` is the DECLARED vertex-axis width: every part is right-padded to it "
    "by repeating its final vertex, and the reader raises rather than silently "
    "widening, so the caller computes it from the frame it is about to write. One row "
    "per shapefile PART, with the record's .dbf attributes replicated."
)

URL_PLACEHOLDER = "REPLACED-BY-THE-CALLER"


# --------------------------------------------------------------------------- io


def load(name: str) -> dict:
    with open(REPO / name) as fh:
        return json.load(fh)


def dump(doc: dict) -> str:
    return json.dumps(doc, indent=1, ensure_ascii=False) + "\n"


# ------------------------------------------------------------------ components


def emission_source(nvert_max: int, geometry: str, stack: bool = False) -> dict:
    """The one uploaded-shapefile loader every generated document reads."""
    return {
        "kind": "static",
        "source": {"url_template": URL_PLACEHOLDER},
        "reference": {
            "notes": (
                f"The caller's own {geometry} emissions, written out of a GeoDataFrame "
                "as a zipped shapefile set and uploaded as an EarthSciLab dataset. The "
                "url is a placeholder: a dispatched run reads the document's own "
                "loaders, so the caller rewrites this with the dataset's store url "
                "before it sends the document. There is no upstream copy to point at, "
                "because the layer is the caller's."
            )
        },
        "metadata": {
            "esio_format": "shapefile",
            "tags": ["shapefile", geometry, "emissions", "user_supplied"],
            "x_esd": {
                "format": "shapefile",
                "container": "zip",
                "crs": "EPSG:4269",
                "note": READER_NOTE,
                # What the caller must put in the .dbf. Stack columns appear
                # only when this document actually reads them.
                "columns": {
                    **{col: "kg/yr" for col in PATHWAYS},
                    **({col: units for col, units, _ in STACK_COLUMNS.values()}
                       if stack else {}),
                },
            },
        },
        "reader_options": {"member": MEMBER, "nvert_max": nvert_max},
        "extent": {"metaparameter": "N_REC"},
    }


def emission_variable(column: str, geometry: str) -> dict:
    """One pollutant's per-record emission, straight off its own column."""
    pathway = PATHWAYS[column]
    return {
        "type": "parameter",
        "units": "kg/yr",
        "default": 0.0,
        "shape": ["emis_records"],
        "description": (
            f"Each {geometry}'s annual {column} emission in kilograms per year, from "
            f"{SOURCE}.{column}. It feeds the {pathway} pathway. Absent from the "
            "caller's frame this column is not written and the pathway is pruned from "
            "the document, so a variable that survives here always has a column behind "
            f"it. `fact` is the kg/yr -> ug/s conversion, applied once to the summed "
            "concentration."
        ),
        "update": {
            "kind": "data",
            "source": SOURCE,
            "from": {"file_variable": column},
        },
    }


def stack_variable(name: str) -> dict:
    column, units, what = STACK_COLUMNS[name]
    return {
        "type": "parameter",
        "units": units,
        "default": 0.0,
        "shape": ["emis_records"],
        "description": (
            f"{what} From {SOURCE}.{column}, already in {units} — the caller's frame "
            "carries SI, so the ft/degF conversions this document's parent needed are "
            "gone and ASME plume rise consumes the column directly."
        ),
        "update": {
            "kind": "data",
            "source": SOURCE,
            "from": {"file_variable": column},
        },
    }


def sr_variable(donor: dict, pathway: str, layer: int) -> dict:
    """An SR matrix for one (pathway, layer), modelled on the donor's own."""
    v = copy.deepcopy(donor)
    v["update"]["source"] = f"ISRM_SR_L{layer}"
    v["update"]["from"]["file_variable"] = pathway
    v["description"] = (
        f"Source-receptor matrix for the {pathway} pathway from emission layer "
        f"{layer}. Rows are the emission-bearing source cells the projection-pushdown "
        f"rewrite derives, columns the N_RCV receptors. Wired from ISRM_SR_L{layer}."
        f"{pathway}; the loader's own gated_select fixes the emission-layer axis, so "
        "this stays rank-2."
    )
    return v


def rename_refs(node, mapping: dict):
    """Rewrite every bare variable reference in an expression tree."""
    if isinstance(node, str):
        return mapping.get(node, node)
    if isinstance(node, list):
        return [rename_refs(x, mapping) for x in node]
    if isinstance(node, dict):
        return {k: rename_refs(v, mapping) for k, v in node.items()}
    return node


def equation(lhs: str, rhs: dict) -> dict:
    """An equation object: `lhs` and `rhs`, and nothing else.

    The schema's equation is a CLOSED object — a `description` key on one is
    rejected outright by `additionalProperties`, which is why every explanation
    in this script lives on the VARIABLE the equation defines.
    """
    return {"lhs": lhs, "rhs": rhs}


def describe(variables: dict, name: str, text: str) -> None:
    """Put prose where the schema allows it: on the variable."""
    variables[name]["description"] = text


def find(equations: list, lhs: str) -> dict:
    for e in equations:
        if e["lhs"] == lhs:
            return e
    raise KeyError(f"no equation for {lhs!r}; have {[e['lhs'] for e in equations]}")


def drop(container, names):
    """Remove `names` from a variables dict or an equations list."""
    names = set(names)
    if isinstance(container, dict):
        for n in names:
            container.pop(n, None)
        return container
    return [e for e in container if e["lhs"] not in names]


def equation_refs(doc: dict) -> dict:
    """`lhs` -> every variable its right-hand side needs, per earthsci-ast.

    `earthsci_ast.free_variables` is used rather than a walk of our own for one
    specific reason: it sees THROUGH `apply_expression_template`. `lcc_forward_x`
    declares `params: ["lon", "lat"]` and its body then reads `lat_1`, `lon_0`,
    `lcc_R`, `lcc_d2r` and `lcc_qp` straight out of the model scope, and
    `krewski_deaths` reads `pop_scale` and `mort_scale` the same way. A walk that
    stops at the call site sees the bindings and none of those, concludes the
    projection constants are unreachable, and prunes them — which passes both the
    schema and `POST /quote`, then projects every source to garbage and dies in
    the fetch with no message at all. Asking the library is how that stays fixed.
    """
    loaded = earthsci_ast.load_document(copy.deepcopy(doc), base_path=str(REPO))
    model = loaded.models[next(iter(doc["models"]))]
    return {eq.lhs: earthsci_ast.free_variables(eq.rhs) for eq in model.equations}


def prune_unreachable(doc: dict, roots) -> list:
    """Keep only what the reported observeds actually need.

    Written as a reachability walk rather than a list of names to delete because
    the alternative does not survive contact with this model: dropping ASME plume
    rise means dropping `stack_layer`, and `stack_layer` is the head of a chain
    through `wind_speed`, `s_class`, `buoy_flux`, `delta_h` and `plume_height`
    that a hand-maintained list gets wrong. A data source no surviving variable
    reads is dropped too — leaving one declared would have the engine discover
    the extent of a file this document no longer looks at.
    """
    model = doc["models"]["ISRM"]
    variables, equations = model["variables"], model["equations"]
    refs = equation_refs(doc)

    needed, stack = set(), [r for r in roots]
    while stack:
        name = stack.pop()
        if name in needed or name not in variables:
            continue
        needed.add(name)
        stack.extend(refs.get(name, set()) & set(variables))

    removed = sorted(set(variables) - needed)
    model["variables"] = {n: v for n, v in variables.items() if n in needed}
    model["equations"] = [e for e in equations if e["lhs"] in needed]

    live = {(v.get("update") or {}).get("source") for v in model["variables"].values()}
    for source in [s for s in doc["data_sources"] if s not in live]:
        doc["data_sources"].pop(source)
        removed.append(f"data_source:{source}")
    return removed


# ------------------------------------------------------------------- pathways


def pathway_layer(doc: dict, donor: dict, columns, layers=(0,), conc_template=None):
    """Give `doc` one conc_/pm_ chain per pollutant, and a summed TotalPM25.

    `donor` is isrm_point.esm's model: the only document that already carries all
    five pathways, so its variables and equations are the source of truth for
    everything except the binning equation, which stays the caller document's.

    `conc_template` is the single-layer contraction, and it has to come from the
    AREA parent rather than from the donor: in the point document
    `conc_PrimaryPM25` is the sum over three emission layers, and only
    `conc_PrimaryPM25_L0` is the matvec. Passing the wrong one wires four
    pathways to a layer sum that does not exist here.
    """
    model = doc["models"]["ISRM"]
    variables, equations = model["variables"], model["equations"]

    sr_donor = donor["variables"]["SR_PrimaryPM25_L0"]
    conc_donor = conc_template if conc_template is not None else find(
        donor["equations"], "conc_PrimaryPM25_L0")

    for column in columns:
        pathway = PATHWAYS[column]
        for layer in layers:
            variables[f"SR_{pathway}_L{layer}"] = sr_variable(sr_donor, pathway, layer)

        if len(layers) > 1:
            # One contraction per emission layer, then the layer sum: exactly what
            # isrm_point.esm does, because plume rise puts mass in all three.
            for layer in layers:
                variables[f"conc_{pathway}_L{layer}"] = copy.deepcopy(
                    donor["variables"][f"conc_{pathway}_L{layer}"])
                equations.append(copy.deepcopy(
                    find(donor["equations"], f"conc_{pathway}_L{layer}")))
            variables[f"conc_{pathway}"] = copy.deepcopy(
                donor["variables"][f"conc_{pathway}"])
            equations.append(copy.deepcopy(find(donor["equations"], f"conc_{pathway}")))
        else:
            variables[f"conc_{pathway}"] = {
                "type": "unknown",
                "units": "ug/m^3",
                "shape": ["rcv_cells"],
                "description": (
                    f"Per-receptor concentration response from the {pathway} pathway: "
                    "the layer-0 SR matrix contracted with the per-cell emissions. One "
                    "contraction rather than the point document's three, because an "
                    "area or line source reaches only emission layer 0."
                ),
            }
            # Both spellings of the donor's emission term are remapped, so the
            # template can come from an area parent (`E_PM25`) or from the point
            # document's layer-0 contraction (`E_PM25_L0`).
            equations.append(equation(
                f"conc_{pathway}",
                rename_refs(copy.deepcopy(conc_donor["rhs"]),
                            {"SR_PrimaryPM25_L0": f"SR_{pathway}_L0",
                             "E_PM25_L0": f"E_{column}",
                             "E_PM25": f"E_{column}"})))

        variables[f"pm_{pathway}"] = {
            "type": "unknown",
            "units": "ug/m^3",
            "shape": ["rcv_cells"],
            "description": f"Reported {pathway} PM2.5 component at each receptor: "
                           f"fact*conc_{pathway}.",
        }
        equations.append(equation(
            f"pm_{pathway}",
            rename_refs(copy.deepcopy(find(donor["equations"], "pm_PrimaryPM25")["rhs"]),
                        {"conc_PrimaryPM25": f"conc_{pathway}"})))

    # TotalPM25 sums whatever pathways survived pruning. Built from the donor's
    # five-term form with the absent terms removed, so the shape of the expression
    # is the point document's and not this script's invention.
    total = copy.deepcopy(find(donor["equations"], "TotalPM25"))
    kept = [PATHWAYS[c] for c in columns]
    total["rhs"]["args"] = [f"conc_{p}" for p in kept]
    total["rhs"]["expr"]["args"][1]["args"] = [
        {"op": "index", "args": [f"conc_{p}", "rcv"]} for p in kept]
    if len(kept) == 1:
        # `+` over one argument is pointless; collapse to the term itself.
        total["rhs"]["expr"]["args"][1] = total["rhs"]["expr"]["args"][1]["args"][0]
    total.pop("description", None)
    total.pop("reference", None)
    describe(variables, "TotalPM25",
             "Total PM2.5 concentration at each receptor: fact times the sum of every "
             f"pathway this document wires ({', '.join(kept)}). This is the quantity "
             "the health functions consume.")
    equations[:] = [e for e in equations if e["lhs"] != "TotalPM25"]
    equations.append(total)
    return doc


def report_block(columns, layers, record_field: str, record_note: str) -> dict:
    """`metadata.x_esd.report`: what a run reports, in the document's vocabulary."""
    pathways = []
    for column in columns:
        pathway = PATHWAYS[column]
        pathways.append({
            "sr_array": pathway,
            "emissions": [f"E_{column}_L{l}" for l in layers] if len(layers) > 1
                         else [f"E_{column}"],
            "concentration": f"conc_{pathway}",
        })
    return {
        "description": (
            "What a run reports, in the document's own vocabulary, so a caller drives "
            "any of these documents without a table of any of them. `pathways` lists "
            "one entry per SR array this document contracts with; `total_pm25` and "
            "`deaths` name the summed observeds; `record_field` names an observed over "
            "the record axis, whose length IS the number of emission records the "
            "loader discovered."
        ),
        "pathways": pathways,
        "total_pm25": "TotalPM25",
        "deaths": {"krewski": "deathsK", "lepeule": "deathsL"},
        "record_field": record_field,
        "record_field_description": record_note,
    }


def head(doc: dict, name: str, description: str, tags, columns, layers,
         record_field, record_note):
    doc["metadata"]["name"] = name
    doc["metadata"]["description"] = description
    doc["metadata"]["tags"] = list(tags)
    doc["metadata"]["x_esd"]["report"] = report_block(
        columns, layers, record_field, record_note)
    doc["metadata"]["x_esd"]["generated_by"] = (
        "data/make_gdf_documents.py — DO NOT EDIT BY HAND. Regenerate after changing "
        "isrm_point.esm, isrm_polygon.esm or isrm_line.esm."
    )
    doc["metadata"]["x_esd"]["column_convention"] = {
        "emissions": {col: "kg/yr" for col in PATHWAYS},
        "stack": {col: units for _, (col, units, _) in STACK_COLUMNS.items()},
        "note": (
            "SI throughout, because the caller hands over a GeoDataFrame rather than a "
            "published inventory in the inventory's own units. An emission column the "
            "frame omits is pruned from the document together with its pathway."
        ),
    }
    return doc


# ------------------------------------------------------------------ documents


def build_point(plume: bool) -> dict:
    """isrm_gdf_point.esm — the point document, read off a shapefile.

    `plume` keeps ASME plume rise and its three emission layers; without it a
    source emits into layer 0 alone and the whole stack apparatus goes, which is
    what a frame with no stack columns can support.
    """
    doc = load("isrm_point.esm")
    model = doc["models"]["ISRM"]
    variables, equations = model["variables"], model["equations"]
    columns = list(PATHWAYS)
    layers = (0, 1, 2) if plume else (0,)

    doc["data_sources"].pop("EGU_Emis")
    doc["data_sources"][SOURCE] = emission_source(1, "point", stack=plume)
    if not plume:
        for layer in (1, 2):
            doc["data_sources"].pop(f"ISRM_SR_L{layer}", None)

    # A pristine copy to lift equations off, taken before anything below edits
    # the working document.
    donor = load("isrm_point.esm")["models"]["ISRM"]

    # FF10's long format is gone, and with it every variable that existed to
    # recover a pathway from a pollutant code. The parent's whole emission and
    # pathway layer goes as well — every E_*, SR_*, conc_* and pm_* is rebuilt
    # below, and leaving the originals in place is how you end up with two
    # equations for E_VOC_L0.
    pathway_wiring = [n for n in list(variables)
                      if n.startswith(("E_", "SR_", "conc_", "pm_"))]
    variables = drop(variables, pathway_wiring +
                     ["pollutant", "emis_annual", "emis_lon", "emis_lat",
                      "is_VOC", "is_NOx", "is_NH3", "is_SOx", "is_PM25",
                      "stkhgt", "stkdiam", "stktemp", "stkvel"])
    equations = drop(equations, pathway_wiring +
                     ["is_VOC", "is_NOx", "is_NH3", "is_SOx", "is_PM25",
                      "stk_height", "stk_diam", "stk_temp", "stk_vel"])

    # The point's own coordinate, out of the shapefile's [record, vertex, xy]
    # geometry. One vertex per record, which is what nvert_max 1 declares.
    doc["metaparameters"]["N_EVERT"] = {
        "type": "integer",
        "default": 1,
        "description": "Vertices per emission record: 1, because a point is one "
                       "vertex. The reader's `nvert_max`.",
    }
    doc["metaparameters"]["N_XY"] = {
        "type": "integer",
        "default": 2,
        "description": "Components of a coordinate: longitude at position 1, latitude "
                       "at position 2. The xy axis of the reader's `geometry`.",
    }
    doc["index_sets"]["emis_vertex"] = {"kind": "interval", "size": "N_EVERT"}
    doc["index_sets"]["xy"] = {"kind": "interval", "size": "N_XY"}
    variables["emis_lonlat"] = {
        "type": "parameter",
        "units": "degree",
        "default": 0.0,
        "shape": ["emis_records", "emis_vertex", "xy"],
        "description": (
            "Every source's coordinate as the shapefile stores it: [record, vertex, "
            "component] in the file's own CRS (EPSG:4269 lon/lat). Component 1 is "
            f"longitude and component 2 latitude. From {SOURCE}.geometry. A point "
            "shapefile has one vertex per record, so the vertex axis is length 1 — the "
            "same geometry variable the area and line documents read, at its "
            "degenerate width."
        ),
        "update": {"kind": "data", "source": SOURCE, "from": {"file_variable": "geometry"}},
    }
    for axis, template in (("X", "lcc_forward_x"), ("Y", "lcc_forward_y")):
        eq = find(equations, axis)
        eq["rhs"] = {
            "op": "aggregate",
            "output_idx": ["r"],
            "ranges": {"r": {"from": "emis_records"}},
            "args": ["emis_lonlat"],
            "expr": {
                "op": "apply_expression_template",
                "args": [],
                "name": template,
                "bindings": {
                    "lon": {"op": "index", "args": ["emis_lonlat", "r", 1, 1]},
                    "lat": {"op": "index", "args": ["emis_lonlat", "r", 1, 2]},
                },
            },
        }
        eq.pop("description", None)
        eq.pop("reference", None)
        describe(variables, axis,
                 f"Each source's projected {axis.lower()} in the ISRM's Lambert "
                 "conformal metres, from the shapefile's single vertex.")

    for column in columns:
        variables[f"emis_{column}"] = emission_variable(column, "point")

    if plume:
        for name in STACK_COLUMNS:
            variables[name] = stack_variable(name)
    # Without plume rise the stack columns are simply never added, and the whole
    # chain that consumed them becomes unreachable. prune_unreachable below is
    # what actually removes it.

    # The binning equations: emis_annual[r]*is_<P>[r] becomes emis_<P>[r], and with
    # no plume rise the layer weight goes too.
    for column in columns:
        for layer in layers:
            eq = copy.deepcopy(find(donor["equations"], f"E_{column}_L{layer}"))
            factors = eq["rhs"]["expr"]["args"][1]["args"]
            keep = [f for f in factors
                    if (f.get("args") or [None])[0] not in ("emis_annual", f"is_{column}")]
            keep.insert(0, {"op": "index", "args": [f"emis_{column}", "r"]})
            eq["rhs"]["expr"]["args"][1]["args"] = keep
            eq["rhs"]["args"] = [a for a in eq["rhs"]["args"]
                                 if a not in ("emis_annual", f"is_{column}")]
            eq["rhs"]["args"].append(f"emis_{column}")
            if not plume:
                eq["rhs"]["expr"]["args"][1]["args"] = [
                    f for f in keep if (f.get("args") or [None])[0] != "w_sr0"]
                eq["rhs"]["args"] = [a for a in eq["rhs"]["args"] if a != "w_sr0"]
            lhs = f"E_{column}_L{layer}" if plume else f"E_{column}"
            eq["lhs"] = lhs
            eq.pop("description", None)
            eq.pop("reference", None)
            prose = (
                f"{column} emissions binned into each source cell"
                + (f", emission layer {layer}: the sum over records contained by the "
                   f"cell of the record's {column} times the plume's layer-{layer} "
                   "weight." if plume else
                   ": the sum over records contained by the cell. The gate ifelse MUST "
                   "stay the FIRST ifelse in this body — the projection-pushdown "
                   "rewrite reads it to derive the support set, and a body it cannot "
                   "recognise silently drops out of `applies_to` and the whole SR slab "
                   "is fetched ungated.")
            )
            variables[lhs] = {
                "type": "unknown",
                "units": "kg/yr",
                "shape": ["src_cells"],
                "description": prose,
            }
            equations.append(eq)

    model["variables"], model["equations"] = variables, equations
    pathway_layer(doc, donor, columns, layers)
    variables["fact"]["default"] = FACT_SI
    variables["fact"]["description"] = (
        "kg/yr -> ug/s: 1e9 ug/kg over 3.1536e7 s/yr. The published documents' "
        "28766.639396245562 is this times 907.18474 kg per short ton, so reading SI "
        "changes this constant and nothing else."
    )

    kind = "with ASME plume rise" if plume else "at ground level"
    return head(
        doc,
        "ISRM_gdf_point" + ("" if plume else "_flat"),
        "The InMAP ISRM over a caller-supplied POINT emission inventory "
        f"({kind}), read from an uploaded shapefile with one column per pollutant "
        "in SI units. Generated from isrm_point.esm by data/make_gdf_documents.py.",
        ["isrm", "inmap", "air_quality", "health", "source_receptor", "point",
         "shapefile", "geodataframe"] + (["plume_rise"] if plume else []),
        columns, layers, "X",
        "An observed over the RECORD axis, whose length IS the number of emission "
        "records the loader discovered. `extent` binds that count to N_REC inside the "
        "engine, but no binding exposes the resolved metaparameter on a prepared "
        "model, so a caller that must report it reads the length of this field.",
    )


def build_area(name: str, geometry: str, record_field: str) -> dict:
    """isrm_gdf_polygon.esm / isrm_gdf_line.esm — one pathway becomes five."""
    doc = load(name)
    model = doc["models"]["ISRM"]
    variables, equations = model["variables"], model["equations"]
    columns = list(PATHWAYS)
    donor = load("isrm_point.esm")["models"]["ISRM"]

    parent_source = "Polygon_Emis" if geometry == "polygon" else "Line_Emis"
    nvert = doc["data_sources"][parent_source]["reader_options"]["nvert_max"]
    doc["data_sources"].pop(parent_source)
    doc["data_sources"][SOURCE] = emission_source(nvert, geometry)

    # Every variable that read the parent loader now reads ours.
    for v in variables.values():
        u = v.get("update") or {}
        if u.get("source") == parent_source:
            u["source"] = SOURCE

    parent_emis = "poly_emis" if geometry == "polygon" else "line_emis"
    binned = "E_PM25"
    # The parent's contraction, taken before the drop below removes it. This is
    # the matvec -- unlike the point document's same-named equation, which sums
    # three emission layers.
    conc_template = copy.deepcopy(find(equations, "conc_PrimaryPM25"))

    if geometry == "line":
        # The caller explodes roads into segments and apportions by length, so the
        # in-document group total goes -- with it an aggregate that ranged over
        # emis_records twice.
        variables = drop(variables, ["line_key", "road_len", "seg_emis"])
        equations = drop(equations, ["road_len", "seg_emis"])
        binned_source = "seg_emis"
    else:
        binned_source = parent_emis
    variables = drop(variables, [parent_emis])

    for column in columns:
        variables[f"emis_{column}"] = emission_variable(column, geometry)

    # The parent's binning equation, and then its whole pathway wiring -- taken
    # and removed BEFORE the new E_* go in, because E_PM25 is both a name the
    # parent used and one this loop is about to create.
    template = copy.deepcopy(find(equations, binned))
    stale = [n for n in list(variables)
             if n.startswith(("SR_", "conc_", "pm_")) or n == binned]
    variables = drop(variables, stale)
    equations = drop(equations, stale)

    for column in columns:
        eq = copy.deepcopy(template)
        eq["lhs"] = f"E_{column}"
        eq["rhs"] = rename_refs(eq["rhs"], {binned_source: f"emis_{column}"})
        eq.pop("description", None)
        eq.pop("reference", None)
        prose = (
            f"{column} emissions binned into each source cell: the sum over records of "
            f"the record's {column} times the fraction of the record that lies in the "
            "cell. The gate ifelse MUST stay the FIRST ifelse in this body — the "
            "projection-pushdown rewrite reads it to derive the support set and the SR "
            "gate, and a body it cannot recognise does not fail, it silently drops out "
            "of `applies_to` and the whole SR slab is fetched ungated."
        )
        variables[f"E_{column}"] = {
            "type": "unknown",
            "units": "kg/yr",
            "shape": ["src_cells"],
            "description": prose,
        }
        equations.append(eq)

    model["variables"], model["equations"] = variables, equations
    pathway_layer(doc, donor, columns, layers=(0,), conc_template=conc_template)
    variables["fact"]["default"] = FACT_SI
    variables["fact"]["description"] = (
        "kg/yr -> ug/s: 1e9 ug/kg over 3.1536e7 s/yr. The published documents' "
        "28766.639396245562 is this times 907.18474 kg per short ton, so reading SI "
        "changes this constant and nothing else."
    )

    allocation = ("the fraction of each polygon's area that lies in the cell"
                  if geometry == "polygon" else
                  "the fraction of each segment's length that lies in the cell")
    return head(
        doc,
        f"ISRM_gdf_{geometry}",
        f"The InMAP ISRM over a caller-supplied {geometry.upper()} emission layer, "
        f"read from an uploaded shapefile with one column per pollutant in SI units "
        f"and allocated to the ISRM source cells by {allocation}. All five pathways, "
        f"lifted off isrm_point.esm. Generated from {name} by "
        "data/make_gdf_documents.py.",
        ["isrm", "inmap", "air_quality", "health", "source_receptor", geometry,
         "shapefile", "geodataframe"],
        columns, (0,), record_field,
        "An observed over the RECORD axis, whose length IS the number of emission "
        "records the loader discovered. `extent` binds that count to N_REC inside the "
        "engine, but no binding exposes the resolved metaparameter on a prepared "
        "model, so a caller that must report it reads the length of this field.",
    )


def validation_errors(doc: dict) -> list:
    """Everything earthsci-ast can say against this document before it is written.

    Worth running here rather than discovering it from the API: EarthSciLab
    rejects an invalid document with a bare `is not valid under any of the
    schemas listed in the 'oneOf' keyword` and 60 KB of echoed model, naming
    nothing. `validate` on the same document says
    `models/ISRM/equations/18: 'description' was unexpected`.
    """
    try:
        loaded = earthsci_ast.load_document(copy.deepcopy(doc), base_path=str(REPO))
        result = earthsci_ast.validate(loaded, base_path=str(REPO))
    except Exception as e:                      # a document too broken to load
        return [f"{type(e).__name__}: {e}"]
    return [str(e) for e in
            list(result.schema_errors) + list(result.structural_errors)]


def vertex_width_errors(doc: dict) -> list:
    """`nvert_max` is two declarations, and they have to agree.

    One is what the reader pads the vertex axis to; the other is how wide the
    document says that axis is. Disagreeing is not a validation error and not a
    wrong number — the engine reads a [records, 5, 2] array into a
    [records, 56, 2] parameter and the worker dies with no message at all.
    """
    model = doc["models"]["ISRM"]
    padded = doc["data_sources"][SOURCE]["reader_options"]["nvert_max"]
    geometry = next((name for name, v in model["variables"].items()
                     if ((v.get("update") or {}).get("from") or {})
                     .get("file_variable") == "geometry"), None)
    if geometry is None:
        return [f"no variable reads {SOURCE}.geometry"]
    axis = model["variables"][geometry]["shape"][1]
    declared = doc["metaparameters"][doc["index_sets"][axis]["size"]]["default"]
    if declared != padded:
        return [f"reader_options.nvert_max is {padded} but {axis} is declared "
                f"{declared} wide"]
    return []


def observable_roots(doc: dict) -> list:
    """What the document promises to be able to report, read off its own report.

    Everything else is reachable from these or is dead weight. `pm_*` are named
    explicitly because `TotalPM25` sums the `conc_*`, not the `pm_*`, so the
    per-pathway components are not reachable from the total.
    """
    report = doc["metadata"]["x_esd"]["report"]
    roots = [report["total_pm25"], *report["deaths"].values(),
             report["record_field"],
             "rcv_W", "rcv_S", "rcv_E", "rcv_N", "rcv_cx", "rcv_cy"]
    for pathway in report["pathways"]:
        roots.append(pathway["concentration"])
        roots.extend(pathway["emissions"])
        roots.append("pm_" + pathway["sr_array"])
    return roots


TARGETS = {
    "isrm_gdf_point.esm": lambda: build_point(plume=True),
    "isrm_gdf_point_flat.esm": lambda: build_point(plume=False),
    "isrm_gdf_polygon.esm": lambda: build_area("isrm_polygon.esm", "polygon", "poly_area"),
    "isrm_gdf_line.esm": lambda: build_area("isrm_line.esm", "line", "seg_len"),
}


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="exit non-zero if any generated document is out of date")
    args = ap.parse_args(argv[1:])

    stale = []
    for name, build in TARGETS.items():
        doc = build()
        removed = prune_unreachable(doc, observable_roots(doc))
        text = dump(doc)
        path = REPO / name
        current = path.read_text() if path.is_file() else None
        if args.check:
            if current != text:
                stale.append(name)
                print(f"  STALE {name}")
            else:
                print(f"  ok    {name}")
            continue
        problems = validation_errors(doc) + vertex_width_errors(doc)
        if problems and not skipped:
            print(f"  INVALID {name} — not written:")
            for problem in problems[:6]:
                print("    ", problem)
            stale.append(name)
            continue
        path.write_text(text)
        model = doc["models"]["ISRM"]
        print(f"  wrote {name}  ({len(text):,} B, {len(model['variables'])} variables, "
              f"{len(model['equations'])} equations, "
              f"{len(doc['data_sources'])} data sources"
              + (f", pruned {len(removed)}" if removed else "") + ", valid)")

    if stale:
        what = "out of date" if args.check else "invalid"
        print(f"\n{len(stale)} document(s) {what}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
