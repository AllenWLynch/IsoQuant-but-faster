"""iqcap: golden-fixture capture harness for IsoQuant 4.1.0 (see SCHEMA.md)."""
import os
import sys

SCHEMA_VERSION = 1

_HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
ISOQUANT_DIR = os.environ.get("IQCAP_ISOQUANT_DIR", os.path.join(PROJECT_ROOT, "IsoQuant"))


def ensure_isoquant_on_path():
    if ISOQUANT_DIR not in sys.path:
        sys.path.insert(0, ISOQUANT_DIR)
