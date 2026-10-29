import struct

import numpy as np

from v1lb.io.bhv2 import read_bhv2


def _var(name, typ, dims, payload=b""):
    out = struct.pack("<Q", len(name)) + name.encode() + struct.pack("<Q", len(typ)) + typ.encode()
    out += struct.pack("<Q", len(dims)) + b"".join(struct.pack("<Q", d) for d in dims)
    return out + payload


def test_bhv2_struct_cell_char_numeric_and_function_handle(tmp_path):
    eye = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])               # 3 x 2, written column-major
    trial = (_var("Trial1", "struct", [1, 1], struct.pack("<Q", 3))
             + _var("TrialError", "double", [1, 1], struct.pack("<d", 0.0))
             + _var("Eye", "double", [3, 2], eye.T.astype("<f8").tobytes())
             + _var("Shapes", "cell", [1, 2], _var("", "char", [1, 3], b"abc") + _var("", "char", [0, 0])))
    fh = _var("fwd", "function_handle", [1, 1], _var("", "char", [1, 14], b"fwd_projective"))
    p = tmp_path / "x.bhv2"
    p.write_bytes(trial + fh)
    d = read_bhv2(p)
    assert d["Trial1"]["TrialError"] == 0.0 and d["Trial1"]["Shapes"] == ["abc", ""]
    np.testing.assert_array_equal(d["Trial1"]["Eye"], eye)
    assert d["fwd"] == "fwd_projective"
