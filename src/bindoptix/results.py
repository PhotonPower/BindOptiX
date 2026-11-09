"""Result envelope ``raytatouille-result`` (Raytatouille ADR 0023), minimal codec.

``{"format": "raytatouille-result", "schema_version": …, "type": …, "data": {…}}``; NumPy arrays
as ``{"dtype", "shape", "values"}`` in C order, complex values as ``[re, im]``, non-finite floats
in arrays as the strings ``"NaN"``, ``"Infinity"``, ``"-Infinity"``. The types written here
(``EmitterResponseTable``, ``FieldMap``) are new and not yet registered in Raytatouille (fusion
plan, decision D13).
"""

from __future__ import annotations

import json
import math
from typing import Any

import numpy as np

FORMAT = "raytatouille-result"
SCHEMA_VERSION = "0.1.0"
_NONFINITE = {"NaN": math.nan, "Infinity": math.inf, "-Infinity": -math.inf}


def _encode_float(x: float) -> float | str:
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "Infinity" if x > 0 else "-Infinity"
    return float(x)


def encode_array(a: np.ndarray[Any, Any]) -> dict[str, Any]:
    """Array → ``{"dtype", "shape", "values"}`` (C order)."""
    arr = np.ascontiguousarray(a)
    flat = arr.ravel()
    if np.iscomplexobj(arr):
        values: list[Any] = [[_encode_float(v.real), _encode_float(v.imag)] for v in flat]
    elif np.issubdtype(arr.dtype, np.floating):
        values = [_encode_float(float(v)) for v in flat]
    elif np.issubdtype(arr.dtype, np.integer) or arr.dtype == np.bool_:
        values = flat.tolist()
    else:
        raise TypeError(f"unsupported dtype {arr.dtype}")
    return {"dtype": str(arr.dtype), "shape": list(arr.shape), "values": values}


def decode_array(d: dict[str, Any]) -> np.ndarray[Any, Any]:
    """Inverse of :func:`encode_array`."""
    dtype = np.dtype(d["dtype"])
    shape = tuple(int(n) for n in d["shape"])
    vals = d["values"]

    def num(v: Any) -> float:
        return _NONFINITE[v] if isinstance(v, str) else float(v)

    if np.issubdtype(dtype, np.complexfloating):
        out = np.array([complex(num(v[0]), num(v[1])) for v in vals], dtype=dtype)
    elif np.issubdtype(dtype, np.floating):
        out = np.array([num(v) for v in vals], dtype=dtype)
    else:
        out = np.array(vals, dtype=dtype)
    if out.size != math.prod(shape):
        raise ValueError("number of values does not match the shape")
    return out.reshape(shape)


def envelope(type_name: str, data: dict[str, Any]) -> dict[str, Any]:
    return {"format": FORMAT, "schema_version": SCHEMA_VERSION, "type": type_name, "data": data}


def to_json(type_name: str, data: dict[str, Any]) -> str:
    return json.dumps(envelope(type_name, data), allow_nan=False, indent=1)


def load_json(text: str, type_name: str) -> dict[str, Any]:
    """Parse an envelope, check format, version 0.1.x and type; return ``data``."""

    def reject(c: str) -> None:
        raise ValueError(f"invalid JSON number literal {c}")

    obj = json.loads(text, parse_constant=reject)
    if obj.get("format") != FORMAT:
        raise ValueError("not a raytatouille-result")
    if not str(obj.get("schema_version", "")).startswith("0.1."):
        raise ValueError(f"unsupported schema_version {obj.get('schema_version')}")
    if obj.get("type") != type_name:
        raise ValueError(f"expected type {type_name}, got {obj.get('type')}")
    data: dict[str, Any] = obj["data"]
    return data
