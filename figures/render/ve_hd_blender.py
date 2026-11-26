"""ve_hd_blender.py (24/11/2026, hình bài CAD R1/R2/R3) — dựa trên bo_dung_chung/n14/ve_n14_blender.py (cùng chất liệu Principled + tối mặt
lưng, cùng ánh sáng hai mặt trời + nền trắng 0,55, cùng camera trực giao nhìn theo +y, cùng mặt cắt nửa: mesh lớn bị cắt trong shader ở
y < Y0, bào quan nhỏ giữ nguyên vẹn nếu trọng tâm đỉnh có y >= Y0). Khác: mỗi bộ phận tô theo HD_i (µm) giữa mesh trước/sau chiếu
bằng LUT tuần tự chung (tệp LUT, thang 0..VMAX); màng trong suốt (ALPHA_mang), ER bán trong suốt (ALPHA_ER), còn lại đặc.
Hình học KHÔNG làm trơn (chỉ smooth shading).
Chế độ (biến môi trường):
  MODE=cat   : mặt cắt nửa toàn tế bào (bỏ phần CLIP_AX < Y0 nếu CLIP_SIGN=lt, mặc định y < 0; hoặc > Y0 nếu gt).
  MODE=zoom  : chỉ bào quan nhỏ (LOAI_NHO) có trọng tâm trong hộp tâm ZC="x,y,z", nửa cạnh ZR (µm); thêm mesh TRƯỚC chiếu
               (TRUOC=<npz>) của chúng dưới dạng khung dây xám đậm + bóng mờ xám; camera trực giao tâm ZC, ortho = 2·ZR.
Biến: FIXED_COL (json loại→hex, không tô HD), MESH (npz sau chiếu), HD (json của hd_n2_dv.py), LUT (json), VMAX, OUT, RES, SPP, ELEV, AZIM, ORTHO, Y0, CX/CZ (tâm khung).
blender -b -P ve_hd_blender.py"""
import bpy, json, os, sys, re, math
import numpy as np
from mathutils import Vector
E = os.environ.get
MODE = E("MODE", "cat"); OUTP = E("OUT"); RES = int(E("RES", "1600")); RESY = int(E("RESY", str(RES))); SPP = int(E("SPP", "256"))
Y0 = float(E("Y0", "0.0")); CAX = E("CLIP_AX", "y"); CSIGN = E("CLIP_SIGN", "lt"); VMAX = float(E("VMAX")); DEV = E("DEV", "GPU")
LUT = json.load(open(E("LUT")))["lut"]; HD = json.load(open(E("HD")))["bo_phan"]
bpy.ops.wm.read_factory_settings(use_empty=True)
sc = bpy.context.scene; sc.render.engine = "CYCLES"
if DEV == "GPU":
    pr = bpy.context.preferences.addons["cycles"].preferences
    for dev in ("OPTIX", "CUDA"):
        try:
            pr.compute_device_type = dev; pr.get_devices()
            if any(d.type == dev for d in pr.devices): break
        except Exception: pass
    for d in pr.devices: d.use = d.type != "CPU"
    sc.cycles.device = "GPU"
else:
    sc.cycles.device = "CPU"
sc.cycles.samples = SPP; sc.cycles.use_denoising = True
try: sc.cycles.denoiser = "OPTIX" if DEV == "GPU" else "OPENIMAGEDENOISE"
except Exception: sc.cycles.denoiser = "OPENIMAGEDENOISE"
sc.cycles.max_bounces = 8; sc.cycles.transparent_max_bounces = 32; sc.render.film_transparent = True
sc.view_settings.view_transform = "Standard"; sc.view_settings.look = "None"; sc.view_settings.exposure = 0.0

def lin(h):
    c = [int(h[i:i + 2], 16) / 255 for i in (1, 3, 5)]
    return tuple(x / 12.92 if x <= 0.04045 else ((x + 0.055) / 1.055) ** 2.4 for x in c)
def mau_hd(h):
    if h is None: return lin("#BBBBBB")
    i = int(round(min(max(h / VMAX, 0.0), 1.0) * (len(LUT) - 1))); return lin(LUT[i])

BIG = set(E("LOAI_LON", "mang,nhan,luoi_noi_chat,golgi,vach,khong_bao").split(","))
LOAI_NHO = set(E("LOAI_NHO", "ti_the,lysosome").split(","))
ALPHA = {"mang": float(E("ALPHA_mang", "0.10")), "vach": float(E("ALPHA_vach", "0.10")),
         "luoi_noi_chat": float(E("ALPHA_ER", "0.35")), "khong_bao": float(E("ALPHA_khong_bao", "0.35"))}

def mat(name, clip, alpha=1.0, rough=0.5, fixed=None, emis=False):
    m = bpy.data.materials.new(name); m.use_nodes = True; nt = m.node_tree; N = nt.nodes; Lk = nt.links
    for n in list(N): N.remove(n)
    out = N.new("ShaderNodeOutputMaterial"); b = N.new("ShaderNodeBsdfPrincipled")
    b.inputs["Roughness"].default_value = rough; b.inputs["Specular IOR Level"].default_value = 0.3
    geo = N.new("ShaderNodeNewGeometry"); oi = N.new("ShaderNodeObjectInfo")
    mix = N.new("ShaderNodeMixRGB")
    if fixed is not None:
        mix.inputs["Color1"].default_value = (*fixed, 1); mix.inputs["Color2"].default_value = (*(0.55 * c for c in fixed), 1)
    else:
        dk = N.new("ShaderNodeMixRGB"); dk.blend_type = "MULTIPLY"; dk.inputs["Fac"].default_value = 1.0
        dk.inputs["Color2"].default_value = (0.55, 0.55, 0.55, 1)
        Lk.new(oi.outputs["Color"], mix.inputs["Color1"]); Lk.new(oi.outputs["Color"], dk.inputs["Color1"]); Lk.new(dk.outputs["Color"], mix.inputs["Color2"])
    Lk.new(geo.outputs["Backfacing"], mix.inputs["Fac"])
    if emis:
        em = N.new("ShaderNodeEmission"); Lk.new(mix.outputs["Color"], em.inputs["Color"]); sh = em.outputs["Emission"]
    else:
        Lk.new(mix.outputs["Color"], b.inputs["Base Color"]); sh = b.outputs["BSDF"]
    if alpha < 1:
        tr = N.new("ShaderNodeBsdfTransparent"); mx = N.new("ShaderNodeMixShader"); mx.inputs["Fac"].default_value = alpha
        Lk.new(tr.outputs["BSDF"], mx.inputs[1]); Lk.new(sh, mx.inputs[2]); sh = mx.outputs["Shader"]
    if clip:
        sep = N.new("ShaderNodeSeparateXYZ"); Lk.new(geo.outputs["Position"], sep.inputs["Vector"])
        lt = N.new("ShaderNodeMath"); lt.operation = "LESS_THAN" if CSIGN == "lt" else "GREATER_THAN"; lt.inputs[1].default_value = Y0
        Lk.new(sep.outputs[CAX.upper()], lt.inputs[0])
        tr2 = N.new("ShaderNodeBsdfTransparent"); mx2 = N.new("ShaderNodeMixShader"); Lk.new(lt.outputs[0], mx2.inputs["Fac"])
        Lk.new(sh, mx2.inputs[1]); Lk.new(tr2.outputs["BSDF"], mx2.inputs[2]); sh = mx2.outputs["Shader"]
    Lk.new(sh, out.inputs["Surface"]); return m

MATS = {}
def mat_for(g):
    key = (g in BIG and MODE == "cat", ALPHA.get(g, 1.0))
    if key not in MATS: MATS[key] = mat("m%d" % len(MATS), key[0], key[1])
    return MATS[key]

def add_mesh(name, V, F):
    me = bpy.data.meshes.new(name); me.vertices.add(len(V)); me.vertices.foreach_set("co", V.ravel())
    me.loops.add(F.size); me.loops.foreach_set("vertex_index", F.ravel())
    me.polygons.add(len(F)); me.polygons.foreach_set("loop_start", np.arange(0, F.size, 3, dtype=np.int32))
    me.update(); me.validate(); me.polygons.foreach_set("use_smooth", np.ones(len(F), bool))
    ob = bpy.data.objects.new(name, me); sc.collection.objects.link(ob); return ob

z = np.load(E("MESH"))
names = sorted(f[:-3] for f in z.files if f.endswith("__v"))
LOAI_MAP = json.loads(E("LOAI_JSON", "{}"))
def loai(n): return LOAI_MAP.get(n) or re.sub(r"_?\d+$", "", n)
BO = set(filter(None, E("BO_LOAI", "").split(",")))
FIX = json.loads(E("FIXED_COL", "{}"))
kept = {}; shown = []
if MODE == "zoom":
    ZC = np.array([float(x) for x in E("ZC").split(",")]); ZR = float(E("ZR"))
for n in names:
    g = loai(n)
    if g in BO: continue
    V = z[n + "__v"].astype(np.float32); F = z[n + "__f"].astype(np.int32); c = V.astype(np.float64).mean(0)
    if MODE == "cat":
        cc = float(c["xyz".index(CAX)])
        if g not in BIG and (cc < Y0 if CSIGN == "lt" else cc > Y0): continue
    else:
        if g not in LOAI_NHO or np.any(np.abs(c - ZC) > ZR): continue
    ob = add_mesh(n, V, F); h = (HD.get(n) or {}).get("hd_max")
    ob.color = (*(lin(FIX[g]) if g in FIX else mau_hd(h)), 1.0); ob.data.materials.append(mat_for(g)); kept[g] = kept.get(g, 0) + 1; shown.append(n)
if MODE == "zoom":
    zt = np.load(E("TRUOC")); wire = mat("wire", False, 1.0, fixed=lin(E("WIRE_COL", "#00A5D9")), emis=True)
    for n in shown:
        if n + "__v" not in zt.files: continue
        V = zt[n + "__v"].astype(np.float32); F = zt[n + "__f"].astype(np.int32)
        ow = add_mesh(n + "_khung", V, F); ow.data.materials.append(wire); ow.visible_shadow = False
        md = ow.modifiers.new("w", "WIREFRAME"); md.thickness = float(E("WIRE_T", "0.008")); md.use_replace = True
        md.use_even_offset = False; md.use_relative_offset = False
print("KEPT", kept, "N", len(shown), flush=True)
w = bpy.data.worlds.new("w"); sc.world = w; w.use_nodes = True
bg = w.node_tree.nodes["Background"]; bg.inputs["Color"].default_value = (1, 1, 1, 1); bg.inputs["Strength"].default_value = 0.55
def sun(rot, e):
    l = bpy.data.lights.new("s", "SUN"); l.energy = e; l.angle = math.radians(12)
    o = bpy.data.objects.new("s", l); sc.collection.objects.link(o); o.rotation_euler = rot
sun((math.radians(55), 0, math.radians(-30)), 2.2)
sun((math.radians(100), 0, math.radians(40)), 0.6)
cd = bpy.data.cameras.new("c"); cd.type = "ORTHO"
cam = bpy.data.objects.new("c", cd); sc.collection.objects.link(cam); sc.camera = cam
EL = math.radians(float(E("ELEV", "15"))); AZ = math.radians(float(E("AZIM", "0")))
if MODE == "cat":
    tgt = Vector((float(E("CX", "0")), 0, float(E("CZ", "0")))); cd.ortho_scale = float(E("ORTHO", "22"))
else:
    tgt = Vector(tuple(ZC)); cd.ortho_scale = float(E("ORTHO", str(2.3 * ZR)))
d = Vector((math.sin(AZ) * math.cos(EL), -math.cos(AZ) * math.cos(EL), math.sin(EL)))
pos = tgt + 60 * d; cam.location = pos
cam.rotation_euler = (tgt - pos).to_track_quat("-Z", "Y").to_euler(); cd.clip_end = 200
sc.render.resolution_x = RES; sc.render.resolution_y = RESY; sc.render.filepath = OUTP
sc.render.image_settings.file_format = "PNG"; sc.render.image_settings.color_mode = "RGBA"
if E("BLEND"): bpy.ops.wm.save_as_mainfile(filepath=E("BLEND"))
bpy.ops.render.render(write_still=True)
print("XONG_VE", OUTP, flush=True)
