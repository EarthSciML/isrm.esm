#!/usr/bin/env python3
"""Validate an ISRM .esm against the EarthSciAST 0.8.0 reference Python binding
and (as a fallback / cross-check) the JSON Schema draft 2020-12.

Usage:
  validate_isrm.py [model.esm ...]     # default: all of isrm*.esm in this dir

The EarthSciAST checkout supplying esm-schema.json is found via $EARTHSCIAST_DIR,
else the sibling ../EarthSciAST checkout. Run under an environment that has the
earthsci_ast package importable (see run-model-py/ for the venv).
"""
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def _earthsciast_dir():
    """Locate the EarthSciAST checkout that ships esm-schema.json."""
    candidates = [
        os.environ.get("EARTHSCIAST_DIR"),
        os.path.join(os.path.dirname(HERE), "EarthSciAST"),
    ]
    for cand in candidates:
        if cand and os.path.isfile(os.path.join(cand, "esm-schema.json")):
            return cand
    raise FileNotFoundError(
        "esm-schema.json not found; set EARTHSCIAST_DIR to an EarthSciAST checkout"
    )


SCHEMA = os.path.join(_earthsciast_dir(), "esm-schema.json")


def reference_binding(ESM):
    import earthsci_ast
    print("== EarthSciAST reference binding ==")
    print("earthsci_ast:", earthsci_ast.__file__)
    esm_file = earthsci_ast.load(ESM)
    print("load: OK")
    result = earthsci_ast.validate(esm_file, base_path=HERE)
    print("is_valid:", result.is_valid)
    print("schema_errors:", len(result.schema_errors))
    for e in result.schema_errors:
        print("  SCHEMA:", e)
    print("structural_errors:", len(result.structural_errors))
    for e in result.structural_errors:
        print("  STRUCT:", e)
    print("unit_warnings:", len(result.unit_warnings))
    for w in result.unit_warnings[:20]:
        print("  UNIT:", w)
    return result.is_valid


def jsonschema_check(ESM):
    import jsonschema
    print("\n== JSON Schema (draft 2020-12) ==")
    with open(SCHEMA) as fh:
        schema = json.load(fh)
    with open(ESM) as fh:
        doc = json.load(fh)
    validator = jsonschema.Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(doc), key=lambda e: list(e.path))
    if not errors:
        print("PASS: document is schema-valid")
        return True
    for e in errors:
        loc = "/".join(str(p) for p in e.path)
        print(f"  ERROR at /{loc}: {e.message}")
    return False


def validate_one(esm):
    print("=" * 72)
    print(esm)
    print("=" * 72)
    ok = True
    try:
        ok = reference_binding(esm) and ok
    except Exception as exc:  # noqa: BLE001
        print("reference binding raised:", type(exc).__name__, exc)
        ok = False
    try:
        ok = jsonschema_check(esm) and ok
    except Exception as exc:  # noqa: BLE001
        print("jsonschema check raised:", type(exc).__name__, exc)
        ok = False
    print(f"\n{os.path.basename(esm)}:", "PASS" if ok else "FAIL", "\n")
    return ok


if __name__ == "__main__":
    targets = sys.argv[1:] or sorted(glob.glob(os.path.join(HERE, "isrm*.esm")))
    if not targets:
        sys.exit("no .esm files to validate")
    print("schema:", SCHEMA)
    results = {esm: validate_one(esm) for esm in targets}
    for esm, ok in results.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {os.path.basename(esm)}")
    print("\nOVERALL:", "PASS" if all(results.values()) else "FAIL")
    sys.exit(0 if all(results.values()) else 1)
