#!/usr/bin/env python3
# =============================================================================
# run_isrm_demo.py — drive the `isrm_demo` analysis of `isrm_point.esm` through
# the EarthSciLab API (https://earthscilab.com) and draw the plot the analysis
# declares.
#
# THIS IS NOT A FOURTH BINDING. `run-rs`, `run-jl` and `run-py` drive the engine
# in-process on this machine; this drives the SAME `earthsci-ast` core through
# the product — on EarthSciLab's hardware, billed to an EarthSciLab account.
# The answer should therefore be the same number, and the totals printed at the
# end exist so you can check that against the records already in this repo.
#
# FOUR THINGS ABOUT THE REQUEST ARE NOT OBVIOUS, and each is a run that fails
# (or silently means nothing) without it:
#
#   1. `"kind": "evaluate"`, not the default `"simulate"`. `models.ISRM` has no
#      `D(·)` anywhere: `system_kind` is `nonlinear`, `isrm_demo`'s time span is
#      0 -> 0, and the whole answer is the observed graph. Dispatched as a
#      simulation the engine does not merely return something meaningless — it
#      refuses the document outright with `Invalid parameter 'src_E'`, which
#      names the wrong thing entirely.
#
#   2. THE TEMPLATE LIBRARY IS INLINED BEFORE THE DOCUMENT IS SENT.
#      `isrm_point.esm` reaches its shared body by `{"ref": "./isrm_base.esm"}`,
#      and a document that arrives over the wire has no directory to be relative
#      to — the server anchors relative refs at ITS working directory, so the
#      raw document fails on the runner with `template-library file not found:
#      isrm_base.esm`. So we merge the library's `expression_templates` into the
#      importing model's own scope and drop the import key.
#
#      Deliberately NOT a full esm-spec §9.7 template resolve, which would also
#      CLOSE the metaparameters — binding `N_REC` to its declared default of 0
#      before the FF10 loader has discovered how many records there are. That
#      run dies inside the engine on a zero-length axis (`col_major_to_arrayd
#      shape mismatch: shape [0] (product 0) vs 1 elements`). `N_REC` is the
#      document's own G6: a metaparameter only the loader knows the value of, so
#      the declaration has to reach the engine still open.
#
#   3. `--records N` IS A DOCUMENT EDIT, not a request parameter. There is no
#      engine-side scale knob to reach for, and that is the right answer — a
#      run's scale is a property of the request, and our requests carry
#      documents. It writes a loader-level `select` range (esm-spec §8.9.2) on
#      every data source that discovers its own extent; because the selection
#      follows the loader's own `record_filter`, `extent` then re-discovers the
#      smaller `N_REC` by itself and nothing else in the document changes.
#
#   4. AUTH IS THE OAuth 2.0 DEVICE GRANT (WorkOS "CLI Auth"), not a pasted
#      token. A full-scale run takes ~50 minutes and a WorkOS access token is
#      short-lived, so anything hand-pasted expires long before the answer
#      exists. The device grant hands back a REFRESH token; this script stores
#      it (0600, `~/.earthscilab/credentials.json`), rotates it, and mints a
#      fresh access token before each request. The token it produces is an
#      ordinary AuthKit user JWT — same JWKS, same `sub` — so the API needed no
#      change to accept it and the run bills to your own account.
#
# WHAT IT PRINTS AND DRAWS
#
#   * the quote (price, machine, cap) and a confirmation prompt;
#   * live progress, with the caveat that the fraction is PHASE-weighted;
#   * sum(TotalPM25), sum(deathsK), sum(deathsL), compared against whichever
#     `run-*/results*.json` in this repo was recorded at the same record count;
#   * the analysis's one declared plot, as a PNG.
#
# USAGE
#
#   python3 run_isrm_demo.py                  # FULL SCALE: ~50 min, ~$0.04
#   python3 run_isrm_demo.py --records 200    # reduced, minutes and pennies
#   python3 run_isrm_demo.py --quote-only     # price it, start nothing (no auth)
#   python3 run_isrm_demo.py --run-id <uuid>  # re-attach to a run already going
#
# Only `matplotlib` is needed beyond the standard library, and only to draw.
# =============================================================================

from __future__ import annotations

import argparse
import base64
import glob
import json
import os
import ssl
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import webbrowser

# --- where things are --------------------------------------------------------

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

API = os.environ.get("EARTHSCILAB_API", "https://api.earthscilab.com")
WORKOS = os.environ.get("WORKOS_API", "https://api.workos.com")
CREDENTIALS = os.environ.get(
    "EARTHSCILAB_CREDENTIALS", os.path.expanduser("~/.earthscilab/credentials.json")
)

MODEL = os.environ.get("ISRM_MODEL", os.path.join(REPO, "isrm_point.esm"))
ANALYSIS = "isrm_demo"

# Fallback only. The names of the things a run reports live in the DOCUMENT, in
# `metadata.x_esd.report` — which is what lets the three shims in this repo name
# no pollutant and no observed of their own — so `totals_for` reads them from
# there, and these are what it uses if a document declares no report block.
#
# They are asked for BY NAME either way, which is also what makes a missing one
# a failed run rather than a quietly short answer: the "give me everything" mode
# skips an observed it cannot evaluate, and this must not.
TOTALS_FALLBACK = ("TotalPM25", "deathsK", "deathsL")

# `GET /datasets/{id}/field` clamps to this and a query cannot raise it. The
# receptor axis is 52,411, so a full field fits and comes back at stride 1; the
# code below checks rather than assuming, because a silently decimated field is
# the one failure this endpoint must not have.
MAX_VALUES = 262_144

# --- palette -----------------------------------------------------------------
# One series, so no legend: the title names it. Blue `#2a78d6` on the light
# chart surface passes the lightness band, the chroma floor and 3:1 contrast.

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
SERIES = "#2a78d6"


# =============================================================================
# HTTP
# =============================================================================


class HttpError(Exception):
    def __init__(self, status: int, body: str, url: str):
        super().__init__(f"HTTP {status} from {url}: {body[:600]}")
        self.status = status
        self.body = body


class OAuthError(Exception):
    """A 400 from WorkOS carrying an OAuth 2.0 `error` code."""

    def __init__(self, error: str, description: str):
        super().__init__(f"{error}: {description}")
        self.error = error


def _open(req: urllib.request.Request, timeout: float):
    try:
        return urllib.request.urlopen(req, timeout=timeout, context=ssl.create_default_context())
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8", "replace")
        try:
            payload = json.loads(body)
        except ValueError:
            payload = {}
        if e.code == 400 and "error" in payload:
            raise OAuthError(payload["error"], payload.get("error_description", "")) from None
        raise HttpError(e.code, body, req.full_url) from None


def http_json(method: str, url: str, *, json_body=None, form=None, headers=None, timeout=120.0):
    data, hdrs = None, dict(headers or {})
    if json_body is not None:
        data = json.dumps(json_body).encode()
        hdrs["Content-Type"] = "application/json"
    elif form is not None:
        data = urllib.parse.urlencode(form).encode()
        hdrs["Content-Type"] = "application/x-www-form-urlencoded"
    req = urllib.request.Request(url, data=data, headers=hdrs, method=method)
    with _open(req, timeout) as resp:
        raw = resp.read()
    return json.loads(raw) if raw else None


# =============================================================================
# Auth — the device grant, and a refresh token that outlives the run
# =============================================================================


class Session:
    """A signed-in EarthSciLab caller, refreshed on demand.

    The whole reason this is a class and not a header constant: a WorkOS access
    token expires in minutes and a full-scale ISRM run takes about an hour, so
    "the token" is a thing that has to be re-derived, not held. Every request
    goes through `headers()`, which mints a new one whenever the current one is
    within a minute of expiry.
    """

    DEVICE_GRANT = "urn:ietf:params:oauth:grant-type:device_code"

    def __init__(self, api: str):
        self.api = api
        self.client_id = http_json("GET", f"{api}/auth/config")["client_id"]
        self._store = self._load()
        self._access = self._store.get("access_token")

    # -- credential file (0600) ----------------------------------------------

    def _load(self) -> dict:
        try:
            with open(CREDENTIALS) as fh:
                return json.load(fh).get(self.api, {})
        except (OSError, ValueError):
            return {}

    def _save(self) -> None:
        try:
            with open(CREDENTIALS) as fh:
                everything = json.load(fh)
        except (OSError, ValueError):
            everything = {}
        everything[self.api] = self._store
        os.makedirs(os.path.dirname(CREDENTIALS), exist_ok=True)
        fd = os.open(CREDENTIALS, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
        with os.fdopen(fd, "w") as fh:
            json.dump(everything, fh, indent=1)

    # -- tokens ---------------------------------------------------------------

    @staticmethod
    def _expiry(token: str | None) -> float:
        """`exp` out of a JWT, read WITHOUT verifying it.

        Reading a claim to decide when to refresh is not the same act as
        trusting one: the API verifies this token against the JWKS, and a lie
        here can only cost us an unnecessary refresh.
        """
        if not token:
            return 0.0
        try:
            payload = token.split(".")[1]
            payload += "=" * (-len(payload) % 4)
            return float(json.loads(base64.urlsafe_b64decode(payload)).get("exp", 0))
        except Exception:
            return 0.0

    def _adopt(self, response: dict) -> None:
        self._access = response["access_token"]
        self._store["access_token"] = self._access
        # Refresh tokens ROTATE. Persisting the new one is not housekeeping —
        # keep the old one and the next run has to sign in again.
        if response.get("refresh_token"):
            self._store["refresh_token"] = response["refresh_token"]
        self._save()

    def _refresh(self) -> bool:
        token = self._store.get("refresh_token")
        if not token:
            return False
        try:
            self._adopt(
                http_json(
                    "POST",
                    f"{WORKOS}/user_management/authenticate",
                    form={
                        "grant_type": "refresh_token",
                        "refresh_token": token,
                        "client_id": self.client_id,
                    },
                )
            )
            return True
        except OAuthError:
            self._store.pop("refresh_token", None)
            return False

    def login(self) -> None:
        """The OAuth 2.0 device authorization grant, start to finish."""
        start = http_json(
            "POST",
            f"{WORKOS}/user_management/authorize/device",
            form={"client_id": self.client_id},
        )
        print()
        print(f"  Sign in to EarthSciLab. Your code is:  {start['user_code']}")
        print(f"  Open: {start['verification_uri_complete']}")
        print()
        try:
            webbrowser.open(start["verification_uri_complete"])
        except Exception:
            pass

        interval = float(start.get("interval", 5))
        deadline = time.time() + float(start.get("expires_in", 300))
        while time.time() < deadline:
            time.sleep(interval)
            try:
                self._adopt(
                    http_json(
                        "POST",
                        f"{WORKOS}/user_management/authenticate",
                        form={
                            "grant_type": self.DEVICE_GRANT,
                            "device_code": start["device_code"],
                            "client_id": self.client_id,
                        },
                    )
                )
                who = http_json("GET", f"{self.api}/me", headers=self.headers())
                print(f"  Signed in as {who.get('email', '?')}.\n")
                return
            except OAuthError as e:
                if e.error == "authorization_pending":
                    continue
                if e.error == "slow_down":
                    interval += 1
                    continue
                raise SystemExit(f"sign-in refused: {e}")
        raise SystemExit("sign-in timed out; run it again")

    def headers(self) -> dict:
        if time.time() > self._expiry(self._access) - 60:
            if not self._refresh():
                self.login()
        return {"Authorization": f"Bearer {self._access}"}

    # -- calls ----------------------------------------------------------------

    def get(self, path: str, timeout: float = 120.0):
        return http_json("GET", f"{self.api}{path}", headers=self.headers(), timeout=timeout)

    def post(self, path: str, body, timeout: float = 300.0):
        return http_json(
            "POST", f"{self.api}{path}", json_body=body, headers=self.headers(), timeout=timeout
        )

    def stream(self, path: str, timeout: float):
        headers = dict(self.headers(), Accept="text/event-stream")
        req = urllib.request.Request(f"{self.api}{path}", headers=headers, method="GET")
        return _open(req, timeout)


# =============================================================================
# The document
# =============================================================================


def inline_template_library(doc: dict, base: str) -> list[str]:
    """Merge every `expression_template_imports` target into its importer.

    See note 2 in the header. The edge in this document is a bare
    `{"ref": "./isrm_base.esm"}` — no bindings, no renames — and the library is
    `esm` + `metadata` + `expression_templates` and nothing else, so a template
    merge is the whole of it. A library carrying index sets, metaparameters or
    per-edge bindings would need more than this.
    """
    inlined = []
    for name, model in (doc.get("models") or {}).items():
        refs = [i["ref"] for i in model.get("expression_template_imports", []) if "ref" in i]
        if not refs:
            continue
        merged: dict = {}
        for ref in refs:
            path = os.path.normpath(os.path.join(base, ref))
            with open(path) as fh:
                merged.update(json.load(fh).get("expression_templates") or {})
            inlined.append(os.path.basename(path))
        model.pop("expression_template_imports", None)
        # The importing document's own declarations win over an imported one.
        merged.update(model.get("expression_templates") or {})
        model["expression_templates"] = merged
    return inlined


def truncate_records(doc: dict, n: int) -> list[str]:
    """Keep only the first `n` records of every self-measuring loader."""
    touched = []
    for name, source in (doc.get("data_sources") or {}).items():
        if (source.get("extent") or {}).get("metaparameter"):
            source["select"] = {"axes": [{"range": {"start": 0, "stop": n}}]}
            touched.append(name)
    return touched


def find_analysis(doc: dict, analysis_id: str) -> tuple[str, dict]:
    for model_name, model in (doc.get("models") or {}).items():
        for analysis in model.get("analyses") or []:
            if analysis.get("id") == analysis_id:
                return model_name, analysis
    have = [
        a.get("id")
        for m in (doc.get("models") or {}).values()
        for a in (m.get("analyses") or [])
    ]
    raise SystemExit(f"no analysis {analysis_id!r} in {MODEL}; it declares {have}")


def totals_for(doc: dict) -> dict[str, str]:
    """The sums to report, named by the document rather than by this file.

    `metadata.x_esd.report` exists so a runner can drive this document and its
    geometry siblings without carrying a table of either: `total_pm25` names one
    observed and `deaths` maps each concentration-response function to its own.
    Reading them here is the same contract `run-py` keeps, one tier out. The key
    is the contract record's own key path, so `oracle` can look the number up
    without a second table.

    `report.record_field` is deliberately NOT requested alongside them. It would
    say how many emission records the loader actually kept — but it lives on the
    record axis rather than the receptor axis, so asking for it makes the run
    write a second, differently-shaped array, and that is an untested output
    shape to put in the path of an hour-long job for a number `--records`
    already determines.
    """
    report = ((doc.get("metadata") or {}).get("x_esd") or {}).get("report") or {}
    named = {}
    if report.get("total_pm25"):
        named["total_pm25"] = report["total_pm25"]
    for function, observed in (report.get("deaths") or {}).items():
        named[f"deaths/{function}"] = observed
    return named or {n: n for n in TOTALS_FALLBACK}


def observeds_for(analysis: dict, totals: dict[str, str]) -> list[str]:
    """What the run has to compute: the totals, plus the plots' own variables.

    Read off the analysis rather than hardcoded, so the request cannot drift
    from the plot it is for — if someone re-points `isrm_demo`'s y axis, this
    asks for the variable the document now names.
    """
    wanted = list(dict.fromkeys(totals.values()))
    for plot in analysis.get("plots") or []:
        for axis in ("x", "y"):
            var = (plot.get(axis) or {}).get("variable")
            if var and var not in wanted:
                wanted.append(var)
    return wanted


# =============================================================================
# The run
# =============================================================================

TERMINAL = {"succeeded", "failed", "cancelled", "capped"}


def money(dollars: float | None) -> str:
    if dollars is None:
        return "—"
    return f"${dollars:.4f}" if 0 < abs(dollars) < 0.01 else f"${dollars:.2f}"


def clock(seconds: float) -> str:
    seconds = int(seconds)
    return f"{seconds // 3600}h{seconds % 3600 // 60:02d}m" if seconds >= 3600 else (
        f"{seconds // 60}m{seconds % 60:02d}s" if seconds >= 60 else f"{seconds}s"
    )


def quote(request: dict) -> dict:
    """`POST /quote` has no database and no auth — it is the pre-login preview."""
    routing = http_json("POST", f"{API}/quote", json_body=request, timeout=300.0)
    option = routing.get("dispatchable")
    if not option:
        raise SystemExit(
            "no dispatchable backend can run a job this size — "
            f"the router says: {routing.get('reason', '(no reason given)')}"
        )
    return routing


def describe_quote(routing: dict, records: int | None) -> float:
    option = routing["dispatchable"]
    estimate = option["estimate"]
    sizing = routing.get("sizing") or {}
    price = estimate["price"]

    print(f"  backend      {option['backend']}  (tier {routing['tier']})")
    if sizing:
        print(
            f"  machine      {sizing.get('vcpus', '?')} vCPU / {sizing.get('memory_mb', '?')} MB"
            f" / {sizing.get('ephemeral_gib', '?')} GiB scratch"
        )
    print(f"  predicted    {clock(estimate['resource_seconds'])}")
    print(f"  cap          {clock(estimate['max_resource_seconds'])}  (stopped and refunded past this)")
    print(f"  queue        ~{clock(estimate['expected_queue_seconds'])}")
    print(f"  price        {money(price)}")
    print(f"  why          {routing.get('reason', '')}")

    # Measured, in docs/isrm.md: the full-scale run that reproduced the oracle
    # took 2,981 s on 4 vCPU / 16 GB. A cap below that is the estimator being
    # wrong in the direction that kills the run — worth saying before you agree
    # to it, even though a capped run is refunded rather than charged.
    if records is None and estimate["max_resource_seconds"] < 3000:
        print()
        print(
            "  ! The cap is under the ~2,981 s a full-scale run has actually taken.\n"
            "    If it fires, the run stops and you are refunded — but you get no answer."
        )
    return price


def watch(session: Session, run_id: str) -> dict:
    """Follow a run to a terminal event, surviving a dropped connection.

    The stream replays everything already recorded before it starts streaming,
    so a reconnect sees what it missed — including a terminal event that landed
    while we were disconnected. That is what makes reconnecting sufficient and a
    separate "did I miss it?" query unnecessary on the common path.
    """
    started = time.time()
    said_phase_caveat = False
    on_a_bar = False

    while True:
        try:
            with session.stream(f"/runs/{run_id}/events", timeout=240.0) as resp:
                for raw in resp:
                    line = raw.decode("utf-8", "replace").strip()
                    if not line.startswith("data:"):
                        continue
                    event = json.loads(line[5:]).get("kind") or {}
                    kind = event.get("type")

                    if kind == "queued":
                        position = event.get("queue_position")
                        print("  queued" + (f" at position {position}" if position else ""), flush=True)
                    elif kind == "started":
                        print("  started", flush=True)
                    elif kind == "progress":
                        fraction = event.get("fraction", 0.0)
                        if fraction > 0.7 and not said_phase_caveat:
                            said_phase_caveat = True
                            print(
                                "\n  (the fraction is PHASE-weighted, and the eight phases are"
                                " nowhere near equal —\n   77% is the S3 fetch 16% in. A long"
                                " crawl here is the I/O, not a hang.)",
                                flush=True,
                            )
                        bar = "#" * int(fraction * 40)
                        print(
                            f"  [{bar:<40}] {fraction * 100:5.1f}%  "
                            f"elapsed {clock(time.time() - started)}",
                            end="\r" if sys.stdout.isatty() else "\n",
                            flush=True,
                        )
                        on_a_bar = sys.stdout.isatty()
                    elif kind in TERMINAL:
                        if on_a_bar:
                            print()
                        return event
        except (HttpError, OAuthError, urllib.error.URLError, OSError, ValueError) as e:
            on_a_bar = False
            print(f"\n  (stream dropped: {e}; the run is server-side and unaffected)", flush=True)

        # The stream ended with no terminal event. Ask outright rather than
        # reporting a failure that did not happen.
        run = session.get(f"/runs/{run_id}")
        if run["status"] in TERMINAL:
            return {"type": run["status"], "message": "(settled by polling)"}
        time.sleep(5)


# =============================================================================
# Reading the answer back
# =============================================================================


def read_fields(session: Session, dataset_id: str, names: list[str]) -> dict[str, list[float]]:
    """One 1-D array per name, out of the run's own Zarr store.

    An evaluate run writes `[eval(1), rcv_cells(52411)]`, so pinning every
    length-1 axis leaves exactly the receptor axis free — "two axes is a field,
    one is a line", and a line is what a scatter wants.
    """
    dataset = session.get(f"/datasets/{dataset_id}")
    pins = ",".join(f"{d['name']}:0" for d in dataset.get("dims", []) if d["size"] == 1)

    series = {}
    for name in names:
        query = {"var": name, "max_values": MAX_VALUES}
        if pins:
            query["at"] = pins
        field = session.get(f"/datasets/{dataset_id}/field?{urllib.parse.urlencode(query)}")
        axes = field["axes"]
        if len(axes) != 1:
            raise SystemExit(
                f"{name} came back with {len(axes)} free axes {[a['name'] for a in axes]}; "
                "expected one. Pin the extra axis and try again."
            )
        if axes[0]["stride"] != 1:
            print(
                f"  ! {name} was DECIMATED to every {axes[0]['stride']}th of "
                f"{axes[0]['stored_size']} — the totals below are of the samples "
                "that came back, not of the field."
            )
        series[name] = field["values"]
    return series


def read_result(session: Session, run_id: str, names: list[str]) -> dict[str, list[float]]:
    """Fallback when the terminal event named no dataset.

    `GET /runs/{id}/result` resolves the run's dataset itself and takes an
    unrestricted slice, flattening a gridded variable to one row per cell keyed
    `name[i]` with 1-BASED indices that are always the store's own — so sorting
    by that index reconstructs the field even if the slice decimated it.
    """
    result = session.get(f"/runs/{run_id}/result", timeout=600.0)
    for entry in result.get("coverage") or []:
        if entry.get("stored_size") and entry.get("size") != entry.get("stored_size"):
            print(f"  ! {entry.get('name')} came back partial: {entry}")

    cells: dict[str, dict[int, float]] = {}
    for row, key in enumerate(result["state_variable_names"]):
        name, _, index = key.partition("[")
        if name in names and index:
            cells.setdefault(name, {})[int(index.rstrip("]").split(",")[0])] = result["state"][row][0]
    return {n: [v for _, v in sorted(cells[n].items())] for n in names if n in cells}


def oracle(
    records: int | None, model_name: str, totals: dict[str, str]
) -> tuple[dict | None, str | None]:
    """The record this repo already has at this record count, if there is one.

    Keyed by `grid.n_rec` and preferring a record of the same document, which is
    how `contract/compare_results.py` matches too. These files are the shims'
    own answers through the same engine, so agreeing with one is evidence the
    product computed the model rather than something adjacent to it.

    `records=None` means the full-scale run, and rather than hardcode what full
    scale IS — the loader discovers `N_REC`, which is the whole point of G6 —
    it takes the largest record count anything here was recorded at.
    """
    candidates = []
    for path in sorted(glob.glob(os.path.join(REPO, "run-*", "results*.json"))):
        try:
            with open(path) as fh:
                record = json.load(fh)
        except (OSError, ValueError):
            continue
        n_rec = (record.get("grid") or {}).get("n_rec")
        # The contract record's key path IS the report key — `total_pm25` and
        # `deaths/<function>` (contract/results_schema.json) — which is why one
        # lookup serves however many response functions a document declares.
        sums = {}
        for key in totals:
            head, _, tail = key.partition("/")
            block = record.get(head) or {}
            value = (block.get(tail) or {}).get("sum") if tail else block.get("sum")
            if value is not None:
                sums[key] = value
        if n_rec is None or not sums:
            continue
        candidates.append((n_rec, record.get("model"), sums, os.path.relpath(path, REPO)))

    if not candidates:
        return None, None
    wanted = records if records is not None else max(n for n, _, _, _ in candidates)
    at_scale = [c for c in candidates if c[0] == wanted]
    if not at_scale:
        return None, None
    # Same document first; otherwise say whose record it is, because the answer
    # moves with the document. `run-rs/results.json` still calls itself
    # `isrm.esm` and docs/isrm.md records that its numbers are the point
    # document's — a mismatch worth printing rather than hiding.
    exact = [c for c in at_scale if c[1] == model_name]
    n_rec, model, sums, path = (exact or at_scale)[0]
    label = f"{path} (n_rec={n_rec:,})"
    if model != model_name:
        label += f" — recorded with model={model!r}, not {model_name!r}"
    return sums, label


# =============================================================================
# The plot
# =============================================================================


def check_can_draw(plots: list[dict]) -> None:
    """Refuse before the run, not after it.

    `draw` imports matplotlib, and importing it for the first time at the END of
    a fifty-minute run would throw the answer away over a missing package. The
    same argument covers the plot type: a form this script cannot draw is worth
    saying now.
    """
    unsupported = [p.get("id") for p in plots if p.get("type") != "scatter"]
    if unsupported:
        raise SystemExit(
            "this script draws scatter only; "
            + ", ".join(repr(i) for i in unsupported)
            + (" are not" if len(unsupported) > 1 else " is not")
        )
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        raise SystemExit(
            "matplotlib is needed to draw the plot, and checking now beats losing a\n"
            "finished run to it. Install it (`pip install -r run-api/requirements.txt`)\n"
            "or pass --quote-only."
        ) from None


def draw(analysis: dict, plot: dict, series: dict[str, list[float]], out: str, subtitle: str):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.transforms import offset_copy

    if plot.get("type") != "scatter":
        raise SystemExit(
            f"plot {plot.get('id')!r} is a {plot.get('type')!r}; this script draws scatter only"
        )

    x_name = plot["x"]["variable"]
    y_name = plot["y"]["variable"]
    x, y = series[x_name], series[y_name]
    if len(x) != len(y):
        raise SystemExit(f"{x_name} has {len(x)} values and {y_name} has {len(y)}")

    fig, ax = plt.subplots(figsize=(9.0, 5.4), dpi=160)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    # 52,411 points on one axis: overplotting, not identity, is the problem to
    # solve. Small marks with low alpha let density read as tone; rasterizing
    # the mark layer alone keeps the file small without softening the text.
    ax.scatter(x, y, s=2.0, c=SERIES, alpha=0.28, linewidths=0.0, rasterized=True)

    ax.set_xlabel(plot["x"].get("label", x_name), color=INK_2, fontsize=9)
    ax.set_ylabel(plot["y"].get("label", y_name), color=INK_2, fontsize=9)
    # Title and subtitle are placed in POINTS above the axes, not in axes
    # fractions: a fraction is a share of the plot's height, so the same 1.02
    # that clears the title on a tall figure sits underneath it on a short one.
    ax.set_title(
        f"{analysis.get('id')} — {y_name} at {len(y):,} ISRM receptor cells",
        color=INK,
        fontsize=12,
        pad=30,
        loc="left",
    )
    ax.text(
        0.0,
        1.0,
        subtitle,
        transform=offset_copy(ax.transAxes, fig=fig, x=0, y=10, units="points"),
        color=MUTED,
        fontsize=8,
        va="bottom",
    )

    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=8, length=3)

    fig.tight_layout()
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)


# =============================================================================
# main
# =============================================================================


def main() -> int:
    parser = argparse.ArgumentParser(
        description=f"Run the {ANALYSIS!r} analysis of {os.path.basename(MODEL)} on EarthSciLab."
    )
    parser.add_argument(
        "--records",
        type=int,
        default=None,
        metavar="N",
        help="keep only the first N emission records (default: FULL SCALE, ~50 min, ~$0.04)",
    )
    parser.add_argument("--esm", default=MODEL, help=f"the document to run (default {MODEL})")
    parser.add_argument("--analysis", default=ANALYSIS, help=f"analysis id (default {ANALYSIS})")
    parser.add_argument("--out", default=os.path.join(HERE, "out"), help="where the PNG goes")
    parser.add_argument("--quote-only", action="store_true", help="price it and stop; no sign-in")
    parser.add_argument("--yes", action="store_true", help="skip the price confirmation")
    parser.add_argument("--run-id", help="re-attach to a run already in flight instead of starting one")
    parser.add_argument("--login", action="store_true", help="sign in and stop")
    parser.add_argument("--logout", action="store_true", help="forget the stored credential and stop")
    args = parser.parse_args()

    if args.logout:
        session_store = {}
        try:
            with open(CREDENTIALS) as fh:
                session_store = json.load(fh)
        except (OSError, ValueError):
            pass
        if session_store.pop(API, None) is not None:
            fd = os.open(CREDENTIALS, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o600)
            with os.fdopen(fd, "w") as fh:
                json.dump(session_store, fh, indent=1)
            print(f"forgot the credential for {API}")
        else:
            print(f"no stored credential for {API}")
        return 0

    if args.login:
        Session(API).login()
        return 0

    # --- the document -------------------------------------------------------

    with open(args.esm) as fh:
        doc = json.load(fh)
    model_name, analysis = find_analysis(doc, args.analysis)
    plots = analysis.get("plots") or []
    if not plots:
        raise SystemExit(f"analysis {args.analysis!r} declares no plots")
    if not args.quote_only:
        check_can_draw(plots)

    inlined = inline_template_library(doc, os.path.dirname(os.path.abspath(args.esm)))
    truncated = truncate_records(doc, args.records) if args.records else []
    totals = totals_for(doc)
    observeds = observeds_for(analysis, totals)

    print(f"\n{os.path.basename(args.esm)} :: {model_name} :: {args.analysis}")
    print(f"  {analysis.get('description', '')}")
    print(f"  inlined      {', '.join(inlined) or '(nothing to inline)'}")
    print(f"  scale        {f'first {args.records:,} records of ' if args.records else 'FULL — all '}"
          f"{', '.join(truncated) if truncated else 'the declared emission inventory'}")
    print(f"  observeds    {', '.join(observeds)}")
    print()

    request = {"esm": doc, "kind": "evaluate", "observeds": observeds}

    # --- price it -----------------------------------------------------------

    price = None
    if not args.run_id:
        print("Quote")
        price = describe_quote(quote(request), args.records)
        print()
        if args.quote_only:
            return 0
        if not args.yes:
            if input(f"Start it for {money(price)}? [y/N] ").strip().lower() not in ("y", "yes"):
                print("nothing started.")
                return 0

    session = Session(API)

    # --- start it, or find it again -----------------------------------------

    if args.run_id:
        run_id = args.run_id
        print(f"Re-attaching to run {run_id}")
    else:
        # Sending back the price that was agreed turns a re-fit between the
        # quote and the dispatch into a refusal instead of a surprise charge.
        run = session.post("/runs", dict(request, max_price=price))
        run_id = run["id"]
        print(f"Run {run_id} — {run['status']} on {run['backend']}, {money(run['price'])}")
        print(f"  re-attach any time with:  --run-id {run_id}")
    print()

    outcome = watch(session, run_id)
    kind = outcome.get("type")
    if kind == "failed":
        raise SystemExit(f"the run failed: {outcome.get('message', '(no message)')}")
    if kind == "cancelled":
        raise SystemExit("the run was cancelled")
    if kind == "capped":
        raise SystemExit(
            f"the run hit its {outcome.get('limit_resource_seconds', 0):.0f}-second cap and was "
            "stopped. Our estimate was too low, so you have not been charged."
        )
    if kind != "succeeded":
        raise SystemExit(f"unexpected terminal event: {outcome}")

    print(f"  succeeded in {clock(outcome.get('resource_seconds', 0))} of resource time")

    # --- read it back -------------------------------------------------------

    dataset_id = outcome.get("dataset_id")
    if dataset_id:
        print(f"  dataset {dataset_id}")
        series = read_fields(session, dataset_id, observeds)
    else:
        print("  (the event named no dataset; falling back to GET /runs/{id}/result)")
        series = read_result(session, run_id, observeds)

    missing = [n for n in observeds if n not in series]
    if missing:
        raise SystemExit(f"the run wrote no {', '.join(missing)}")

    # --- the totals ---------------------------------------------------------

    print("\nTotals")
    reference, source = oracle(args.records, os.path.basename(args.esm), totals)
    for key, name in totals.items():
        total = sum(series[name])
        line = f"  sum({name})".ljust(20) + repr(total)
        expected = (reference or {}).get(key)
        if expected:
            rel = abs(total - expected) / abs(expected)
            line += f"   vs {expected!r}   rel {rel:.1e}"
        print(line)
    if source:
        print(f"  reference: {source}")
    else:
        print("  (no record in this repo at this record count to compare against)")

    # --- the plot -----------------------------------------------------------

    os.makedirs(args.out, exist_ok=True)
    scale = f"first {args.records:,} records" if args.records else "full inventory"
    print()
    for plot in plots:
        path = os.path.join(args.out, f"{args.analysis}_{plot['id']}.png")
        try:
            draw(
                analysis,
                plot,
                series,
                path,
                f"EarthSciLab run {run_id} · {os.path.basename(args.esm)} · {scale}",
            )
        except Exception as e:
            # The numbers cost real time and real money; a drawing bug must not
            # be what loses them. Spill them and re-raise.
            spill = os.path.join(args.out, f"{args.analysis}_{run_id}.json")
            with open(spill, "w") as fh:
                json.dump(series, fh)
            raise SystemExit(f"drawing {plot['id']} failed ({e}); the values are in {spill}")
        print(f"  wrote {path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
