# -*- coding: utf-8 -*-
"""igl_tho.py (22/11/2026, v123u5) — tiến trình phụ cho bước D_i của botri_cucbo_v123u5: nhận (P float64 (n,3), V float64, F int64) qua stdin,
trả igl.signed_distance(P, V, F)[0] (float64) qua stdout. Cùng hàm, cùng đầu vào như mã gốc; chỉ chạy trên lõi khác (taskset của bên gọi)."""
import sys, numpy as np, igl
i_, o_ = sys.stdin.buffer, sys.stdout.buffer
def doc(n):
    b = bytearray()
    while len(b) < n:
        x = i_.read(n - len(b))
        if not x: return None
        b += x
    return bytes(b)
while True:
    h = doc(24)
    if h is None: break
    n, nv, nf = np.frombuffer(h, np.int64).tolist()
    P = np.frombuffer(doc(24 * n), np.float64).reshape(n, 3); V = np.frombuffer(doc(24 * nv), np.float64).reshape(nv, 3)
    F = np.frombuffer(doc(24 * nf), np.int64).reshape(nf, 3)
    d = igl.signed_distance(P, V, F)[0]
    o_.write(np.ascontiguousarray(d, np.float64).tobytes()); o_.flush()
