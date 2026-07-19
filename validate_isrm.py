#!/usr/bin/env python3
"""Validate isrm.esm against the EarthSciAST 0.8.0 reference Python binding
and (as a fallback / cross-check) the JSON Schema draft 2020-12.

Run with the earthsci-ast-py venv, e.g.:
  /Users/ctessum/code/earthsciml/EarthSciAST/pkg/earthsci-ast-py/.venv/bin/python3 validate_isrm.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ESM = os.path.join(HERE, "isrm.esm")
SCHEMA = "/Users/ctessum/code/earthsciml/EarthSciAST/esm-schema.json"


def reference_binding():
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


def jsonschema_check():
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


if __name__ == "__main__":
    ok = True
    try:
        ok = reference_binding() and ok
    except Exception as exc:  # noqa: BLE001
        print("reference binding raised:", type(exc).__name__, exc)
        ok = False
    try:
        ok = jsonschema_check() and ok
    except Exception as exc:  # noqa: BLE001
        print("jsonschema check raised:", type(exc).__name__, exc)
        ok = False
    print("\nOVERALL:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)
