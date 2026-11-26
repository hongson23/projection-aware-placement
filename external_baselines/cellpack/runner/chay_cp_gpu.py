"""chay_cp_gpu.py — chạy cellPACK trên đầu vào do cau_noi_v12b.py xuat ghi ra, GIỐNG HỆT chay_cellpack_b.py
(cùng recipe, cùng config, cùng thứ tự gọi, cùng bộ đếm tiến độ), chỉ khác: đọc đầu vào từ --vao (chỉ đọc), ghi mọi thứ vào --ra_dir,
tuỳ chọn cProfile (--prof, SIGUSR1 ghi hồ sơ giữa chừng), và ghi thêm dấu vân tay để so tương đương:
  - tu_the_cp.json: như bản gốc (dat {ten: {t, R}}), số thực ghi bằng repr (17 chữ số) => so từng bit;
  - van_tay.json: sha256 của lưới sau buildGrid (vị trí, khoảng cách, id ngăn), của trạng thái lưới sau pack_grid,
    của trạng thái bộ sinh số ngẫu nhiên (numpy + random) sau pack, và của danh sách đặt (thứ tự, tên, vị trí, xoay).
Bản cellPACK nào được nạp do PYTHONPATH quyết (gốc: không đặt; GPU: cellpack_gpu/src). Ghi đường dẫn cellpack.__file__ ra để kiểm.
Dùng: python chay_cp_gpu.py --vao VAO --ra_dir RA --seed S [--thu_tu ban_kinh_giam] [--bo tac_gia] [--prof]"""
import argparse, json, os, sys, time, subprocess, hashlib
ap = argparse.ArgumentParser(); ap.add_argument("--vao", required=True); ap.add_argument("--ra_dir", required=True)
ap.add_argument("--seed", type=int, required=True)
ap.add_argument("--bo", default="tac_gia", choices=["mac_dinh", "tac_gia"]); ap.add_argument("--spacing", type=float, default=0.0)
ap.add_argument("--thu_tu", default="ban_kinh_giam", choices=["ban_kinh_giam", "mac_dinh"])
ap.add_argument("--place_method", default="jitter", choices=["jitter"]); ap.add_argument("--ra", default="tu_the_cp.json")
ap.add_argument("--luoi_tu", default=None, help="CHỈ ĐỂ PHÁT TRIỂN: nạp lưới đã lưu (pickle cellPACK) thay vì dựng lại");
ap.add_argument("--prof", action="store_true"); ap.add_argument("--chi_luoi", action="store_true", help="chỉ dựng lưới rồi dừng (thử tương đương lưới)")
a = ap.parse_args(); V = os.path.abspath(a.vao); D = os.path.abspath(a.ra_dir); os.makedirs(D, exist_ok=True)
if a.prof:
    import cProfile, signal
    _PR = cProfile.Profile()
    def _dump(*_):
        _PR.disable(); _PR.dump_stats(f"{D}/prof.pstats"); _PR.enable(); print("PROF_GHI", time.strftime("%T"), flush=True)
    signal.signal(signal.SIGUSR1, _dump); _PR.enable()
dv = json.load(open(f"{V}/dau_vao_cellpack.json"))
JA = {"mac_dinh": 5, "tac_gia": 300}[a.bo]; RT = {"mac_dinh": 30, "tac_gia": 30}[a.bo]
obj = {}; comp = {}
for k in ("vo_chua", "vung_cam", "long_nhan"):
    obj[k + "_mesh"] = {"type": "mesh", "representations": {"mesh": {"path": V, "name": dv["ngan"][k]["tep"], "format": "obj"}}}
noi = {"te_bao_chat": [], "long_nhan": []}
for b in dv["bo_phan"]:
    obj[b["ten"]] = {"type": "multi_sphere", "jitter_attempts": JA, "rejection_threshold": RT, "place_method": a.place_method, "packing_mode": "random",
                     "rotation_range": 6.2831, "use_rotation_axis": False, "cutoff_surface": 1e-6, "max_jitter": [1, 1, 1],
                     **({"priority": -round(float(b["R0"]), 6)} if a.thu_tu == "ban_kinh_giam" else {}),
                     "representations": {"packing": {"path": f"{D}/sph", "name": b["ten"] + ".sph", "format": "sph"}}}
    os.makedirs(f"{D}/sph", exist_ok=True)
    with open(f"{D}/sph/{b['ten']}.sph", "w") as fs:
        n_ = len(b["r"]); fs.write(f"{min(b['r']):.6f} {b['R0']:.6f}\n2\n1\n0.0 0.0 0.0 {b['R0']:.6f} " + " ".join(map(str, range(n_))) + f"\n{n_}\n")
        for (x, y, z), r in zip(b["tam"], b["r"]): fs.write(f"{x:.6f} {y:.6f} {z:.6f} {r:.6f}\n")
    noi[b["ngan"]].append({"object": b["ten"], "count": 1})
comp["bounding_area"] = {"regions": {"interior": ["te_bao"]}}
comp["te_bao"] = {"object": "vo_chua_mesh", "count": 1, "regions": {"interior": ["vung_cam"] + noi["te_bao_chat"]}}
comp["vung_cam"] = {"object": "vung_cam_mesh", "count": 1, "regions": {"interior": ["long_nhan"]}}
comp["long_nhan"] = {"object": "long_nhan_mesh", "count": 1, "regions": {"interior": noi["long_nhan"]}}
rec = {"name": f"cau_noi_{a.seed}", "version": "1.0", "format_version": "2.1", "bounding_box": dv["hop"], "objects": obj, "composition": comp}
cfg = {"name": "cau_noi", "out": f"{D}/cp_out", "place_method": a.place_method, "overwrite_place_method": False, "save_analyze_result": False,
       "number_of_packings": 1, "show_progress_bar": False, "randomness_seed": a.seed, "clean_grid_cache": True, "load_from_grid_file": False}
if a.spacing > 0: cfg["spacing"] = a.spacing
json.dump(rec, open(f"{D}/recipe_cp.json", "w"), indent=1); json.dump(cfg, open(f"{D}/config_cp.json", "w"), indent=1)
import numpy as np, random
from cellpack import autopack
from cellpack.autopack import upy
from cellpack.autopack.Environment import Environment
from cellpack.autopack.loaders.config_loader import ConfigLoader
from cellpack.autopack.loaders.recipe_loader import RecipeLoader
import logging; logging.getLogger().setLevel(logging.WARNING)
for n_ in ("cellpack", "autopack"): logging.getLogger(n_).setLevel(logging.WARNING)
print("CELLPACK_FILE", autopack.__file__, flush=True)
def _h(*arrs):
    h = hashlib.sha256()
    for x in arrs:
        x = np.ascontiguousarray(np.asarray(x)); h.update(str(x.dtype).encode()); h.update(str(x.shape).encode()); h.update(x.tobytes())
    return h.hexdigest()
def _van_tay_luoi(env):
    g = env.grid; out = {}
    for nm in ("masterGridPositions", "distToClosestSurf", "compartment_ids", "freePointMask", "free_points"):
        if hasattr(g, nm):
            try: out[nm] = _h(getattr(g, nm))
            except Exception as e: out[nm] = f"loi {e}"
    out["nbFreePoints"] = int(getattr(g, "nbFreePoints", -1))
    for c in env.compartments:
        for nm in ("insidePoints", "surfacePoints", "surfacePointsNormals", "surfacePointsCoords"):
            if hasattr(c, nm):
                try: out[f"{c.name}.{nm}"] = _h(getattr(c, nm) if not isinstance(getattr(c, nm), dict) else sorted(getattr(c, nm).items()))
                except Exception as e: out[f"{c.name}.{nm}"] = f"loi {e}"
    return out
cfgd = ConfigLoader(f"{D}/config_cp.json").config; recd = RecipeLoader(f"{D}/recipe_cp.json", False).recipe_data
helper = upy.getHelperClass()(vi="nogui"); autopack.helper = helper
env = Environment(config=cfgd, recipe=recd); env.helper = helper
if a.luoi_tu: env.previous_grid_file = a.luoi_tu
t0 = time.time(); env.buildGrid(rebuild=True); t_luoi = time.time() - t0
vt = {"cellpack_file": autopack.__file__, "luoi": _van_tay_luoi(env), "giay_luoi": t_luoi}
print('LUOI_XONG', round(t_luoi, 1), 's spacing', env.grid.gridSpacing, 'diem', len(env.grid.masterGridPositions), flush=True)
json.dump(vt, open(f"{D}/van_tay.json", "w"), indent=1)
if a.chi_luoi:
    if a.prof: _PR.disable(); _PR.dump_stats(f"{D}/prof.pstats")
    print("CHI_LUOI_XONG", flush=True); sys.exit(0)
from cellpack.autopack.ingredient.Ingredient import Ingredient as _Ing
_DEM = {"dat": 0, "loai": 0, "t": time.time()}
_place0, _reject0 = _Ing.place, _Ing.reject
def _place(self, *x, **y):
    _DEM["dat"] += 1
    if _DEM["dat"] % 25 == 0: print("TIEN_DO dat", _DEM["dat"], "/", len(dv["bo_phan"]), "loai", _DEM["loai"], "giay", round(time.time() - _DEM["t"], 1), flush=True)
    return _place0(self, *x, **y)
def _reject(self, *x, **y):
    _DEM["loai"] += 1; return _reject0(self, *x, **y)
_Ing.place, _Ing.reject = _place, _reject
t1 = time.time(); _DEM["t"] = t1; env.pack_grid(seedNum=a.seed, verbose=0, usePP=False); t_dat = time.time() - t1
dat = {}; thu_tu = []
for o in env.packed_objects.get_ingredients():
    if o.is_compartment: continue
    R = np.asarray(o.rotation, float)[:3, :3]; ten = o.name.split("_interior_", 1)[-1]
    dat[ten] = {"ten_cp": o.name, "t": np.asarray(o.position, float).tolist(), "R": R.tolist()}
    thu_tu.append([ten, np.asarray(o.position, float).tolist(), R.tolist()])
vt["sau_dat"] = _van_tay_luoi(env); vt["giay_dat"] = t_dat
vt["dat_sha"] = hashlib.sha256(json.dumps(thu_tu).encode()).hexdigest()
vt["np_rng_sha"] = hashlib.sha256(repr(np.random.get_state()).encode()).hexdigest()
vt["py_rng_sha"] = hashlib.sha256(repr(random.getstate()).encode()).hexdigest()
vt["so_dat"] = len(dat); vt["so_lan_loai"] = _DEM["loai"]
try:
    from cellpack.autopack import gpu_tang_toc as _gt
    vt["gpu_thong_ke"] = dict(_gt.THONG_KE); vt["gpu_bat"] = _gt.BAT; vt["gpu_kiem"] = _gt.KIEM; vt["gpu_dung"] = _gt.dung_gpu()
    _L = _gt._LUOI["obj"]; vt["gpu_mang_luoi_dung"] = None if _L is None else _L.mang_luoi_dung
    print("GPU_THONG_KE", vt["gpu_thong_ke"], "luoi_dung", vt["gpu_mang_luoi_dung"], flush=True)
except ImportError:
    vt["gpu_thong_ke"] = None
json.dump(vt, open(f"{D}/van_tay.json", "w"), indent=1)
try: rev = subprocess.run(["git", "-C", os.path.join(os.path.dirname(autopack.__file__), "..", ".."), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
except Exception: rev = None
out = {"md5_dau_vao": hashlib.md5(open(f"{V}/dau_vao_cellpack.json", "rb").read()).hexdigest(), "dat": dat, "giay_dat": round(t_dat, 3), "giay_luoi": round(t_luoi, 3), "so_yeu_cau": len(dv["bo_phan"]), "so_dat_cp": len(dat),
       "cellpack_commit": rev, "cellpack_file": autopack.__file__, "bo": a.bo, "thu_tu": a.thu_tu, "jitter_attempts": JA, "rejection_threshold": RT, "place_method": a.place_method, "seed_cp": a.seed,
       "so_lan_loai": _DEM["loai"], "spacing_luoi": float(env.grid.gridSpacing), "so_diem_luoi": int(len(env.grid.masterGridPositions)), "nguong": dv["nguong"]}
json.dump(out, open(f"{D}/{a.ra}", "w"), indent=1)
if a.prof: _PR.disable(); _PR.dump_stats(f"{D}/prof.pstats")
print("CP_XONG dat", len(dat), "/", len(dv["bo_phan"]), "giay_luoi", round(t_luoi, 1), "giay_dat", round(t_dat, 1), "spacing", out["spacing_luoi"], "diem", out["so_diem_luoi"], flush=True)
