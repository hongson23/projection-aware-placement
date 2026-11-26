# PA-PLACE: projection-aware placement for dense 3D packing of nested parts under containment and clearance constraints

Code, frozen rule files, cell specifications and reviewed results for the paper

> H. S. Hoang, T. T. Dang, D. H. Nguyen. *Projection-aware placement for dense 3D packing of nested parts under containment and clearance constraints.* Submitted to Computer-Aided Design, 2026.

The pipeline places the parts of a cell (organelles inside a membrane, a nucleus inside the cell, sub-parts inside parts) as signed distance fields, then runs PROJECT, a projection of all fields onto the rule set (containment, clearance and contact rows), meshes the projected fields, applies one uniform mesh repair to every method, and scores the result with an independent rule checker. PLACE is the gradient-based placement stage. **PA-PLACE** adds a projection-aware term that predicts, during placement, how much PROJECT will carve each part and ranks restarts by the predicted deviation.

The source files are the versions that produced the numbers in the paper, with their comments removed for this release; nothing else was changed. Every Python file has the same abstract syntax tree as the version that was run, the shell scripts differ only by removed comment lines, and the macro file defines the same values (`CODE_EQUIVALENCE.tsv` lists, per file, the md5 of the version that was run and of the released file). `MANIFEST.md5` lists the released files (`md5sum -c MANIFEST.md5`). Identifiers, docstrings and printed messages are in Vietnamese; the glossary below maps the internal names to the names used in the paper.

## Repository layout

| directory | contents |
|---|---|
| `builder/` | the builder: PLACE, the projection-aware term, PROJECT, meshing (`botri_chung_*.py`), the local-grid module used for plant cells (`botri_cucbo_*.py`), and an optional helper process (`igl_tho.py`) |
| `checker/` | independent checker v3 (`cham_mot.py` scores one cell), mesh sanitation (`ve_sinh_luoi.py`), the uniform mesh repair (`sua_mesh_dn.py`), mutation tests (`cay_loi*.py`, `tong_hop_cong_v3b.py`), the rule copies it reads (`luat_v12/`), and the frozen md5 lists (`DONG_BANG_BO_KIEM_V3.txt`, `VE_SINH_LUOI_MD5.txt`, `SUA_MESH_MD5.txt`) |
| `rules/` | frozen fact-rule files v1.2 (animal 20 rules, plant 28 rules, 152 convention parameters) and their SHA-256 freeze record |
| `specs/` | cell specifications: `v12/` synthetic animal and plant cells, `n14/` the five real Allen hiPSC containers, their preparation script and the D5-like organelle counts |
| `analysis/` | scripts that turned per-run outputs into the reviewed result files |
| `external_baselines/` | bridges to cellPACK and Packmol, the GPU patch of cellPACK, and the cellPACK license |
| `figures/` | Blender render scripts for the 3D figures and their description (`HINH_CAD.md`) |
| `results/` | reviewed numbers (JSON, CSV, MD) behind every table, the macro generator `so_lieu/tao_so_lieu.py`, the 2D paired-seed figure script, mutation-test results, and `commands/paper_runs.tsv` with the exact command of each of the 513 runs |
| `environment/` | frozen package lists and hardware description |

## Glossary of internal names

| internal name | meaning in the paper |
|---|---|
| `botri_chung_*.py` (bố trí chung) | builder |
| `bo_kiem` (bộ kiểm), `kiem3.py` | checker (v3) |
| `dac_ta` (đặc tả) | cell specification |
| `luat_su_that` / `tham_so_quy_uoc` | fact rules / convention parameters of the frozen rule file |
| `khuon` (khuôn), `--thu_vien T` / `A` | template; procedural / AI-generated template library |
| `U5`, `--u5` | projection-aware term |
| `U5b` (run ids `*_U5b_*`) | **PA-PLACE** |
| `td` (tự động), `--cach tu_dong` | PLACE |
| `P4`, `P8` | PLACE with K = 4 / K = 8 restarts |
| `U5b8` | PA-PLACE with K = 8 |
| `nn1x` (ngẫu nhiên), `--cach ngau_nhien` | rejection sampling (1x budget) |
| `nn1xm` | rejection sampling keeping the candidate with the best soft objective (Rej.+soft, TV-D1) |
| `un1x`, `unc8` (ủ nhiệt), `--cach u_nhiet` | simulated annealing (SA) |
| `u5ktt`, `--ktt` (khả thi trước) | PA-PLACE with the feasibility-first variant (PA-PLACE+F) |
| `a1`, `l0`, `ah`, `wv`, `a1l03`, `a1l3` | mechanism controls: a ≡ 1, λ = 0, a ≡ ½, W ∝ V, a ≡ 1 with λ = 0.3 / 3 |
| `CXe`, `CX2e`, `CX2b`, `--cx` | ComPLx-style in-loop anchoring baseline CX and its versions |
| `cellpack` / `cp`, `packmol` / `pk` | external baselines |
| `k08` (`dong_vat_v12_khe0.8.json`) | K08: D5 with the inter-kind clearance raised to 0.8 µm |
| `D1`, `D4`, `D5`, `D4A2` | animal configurations |
| `TV-D1`, `TV-D5`, `TV-real` (thực vật) | plant configurations |
| `N14`, `n14d5` | real Allen hiPSC containers at D5-like density |
| `N1` | external-baseline comparison (Table tab:ext) |
| `W1` | main animal run (Table tab:w1) |
| `kq` (kết quả), `cham` (chấm) | run outputs, scoring |
| `mesh_vo_bao.npz`, `mesh_vo_bao_truoc.npz`, `thong_tin.json` | final meshes, pre-projection meshes, run report of one cell |
| `ER` (lưới nội chất), `nhan`, `ti_the`, `lysosome`, `luc_lap`, `golgi`, `peroxisome` | ER, nucleus, mitochondria, lysosomes, chloroplasts, Golgi, peroxisomes |

Main builder flags: `--cach` placement method, `--K` restarts, `--so TYPE=N` organelle count, `--tang2` build tier-2 sub-parts, `--khong_tu_do` leave all mesh measurements to the checker, `--eval_dung_tu_dau` exact far-field evaluation from the start, `--hd_truoc_sau` measure the pre/post-projection Hausdorff distance, `--u5_lambda` λ, `--u5_kappa` hinge κ, `--u5_chon gc` rank restarts by the predicted Ĝ, `--u5_er` charge the dependent region (ER), `--u5_tam` add the deepest-node sample. The help string of every flag is in the `argparse` block at the top of each builder.

## Which builder produced which result

| paper item | builder | example run id (full command in `results/commands/paper_runs.tsv`) | analysis → reviewed result |
|---|---|---|---|
| Sec. 4 (PROJECT, PLACE, projection-aware term, feasibility-first variant) | `builder/botri_chung_v122u5e.py`; plant: `builder/botri_chung_v123u5c.py` + `botri_cucbo_v123u5b.py` | — | — |
| Sec. 4.5 (meshing, uniform repair, checker, mutation tests) | `checker/ve_sinh_luoi.py`, `checker/sua_mesh_dn.py`, `checker/cham_mot.py`; mutation tests `checker/cay_loi.py`, `cay_loi_tv.py`, `v3_cay_loi_chung.py`, `dung_cong_v3.sh`, `tong_hop_cong_v3b.py` | — | `results/checker_mutation_tests/cong_v3_cuoi.json` |
| Sec. 6.1, Table tab:w1 (K08, D5, D4A2; seeds 6511–6520) | `botri_chung_v122u5e.py` | `u5_r_U5b_k08_6511`, `u5_r_P4_k08_6511`, `u5_r_P8_k08_6511` | `analysis/tang_dg/u5v2/phan_tich_doc_lap_u5.py` → `phan_tich_rong.py` → `results/so_lieu/nguon/ket_qua_rong.json` |
| Fig. fig:paired (per-seed pairs, K08) | — | — | `results/hinh/fig4_paired_seeds.py` (reads `ket_qua_cx2.json`) |
| Fig. fig:k08 (3D render, seed 6512) | — | `u5_r_P8_k08_6512`, `u5_r_U5b_k08_6512` | `figures/render/ve_r1.sh`, `cfg_r1.py`, `hd_n2_dv.py`, `ve_hd_blender.py`, `ghep_v2.py`; `figures/HINH_CAD.md` |
| Sec. 6.2, Table tab:abl (mechanism controls) | `botri_chung_v122u5f.py` (a1l03, a1l3, ah, wv); `botri_chung_v122u5e.py` (a1, l0) | `u5_c_a1l03_k08_6511`, `u5_r_l0_k08_6511` | `analysis/tang_dg/u5v2/phan_tich_u5c.py` → `ket_qua_u5c.json` |
| Sec. 6.3, CX comparison (Supp. S2) | `botri_chung_v122u5x.py` | `u5_x_CXe_k08_6511`, `u5_x_CX2e_k08_6511`, `u5_x_CX2b_k08_6511` | `analysis/tang_dg/u5x/phan_tich_cx2.py` (+ `dem_buoc.py`, `chon_dev2.py`) → `ket_qua_cx2.json`, `chon_dev2.json` |
| Sec. 6.4, Table tab:plant, Fig. fig:plant (seeds 7311–7316; TV-D1 7311–7313) | `botri_chung_v123u5c.py` | `u5_tv_c_real_u5ktt_7311`, `u5_tv_c_real_td_7311`, `u5_tv_c_real_nn1x_7311`, `u5_tv_c_real_un1x_7311` | `analysis/tang_dg/u5tv/n2_tv.py` (HD measure), `phan_tich_tv_u5.py` → `ket_qua_tv.json`; render `figures/render/ve_r2.sh` |
| Sec. 6.5, Table tab:ext (D1, D4, D4A2, D5; seeds 6201–6210; Supp. S4) | `botri_chung_v122u5f.py` (PA-PLACE), `botri_chung_v122e3.py` (PLACE, rejection, SA, and PROJECT for cellPACK/Packmol poses) | `u5_c_U5b_D5_6201`, `e3b_B_D5_td_6201`, `e3b_B_D5_nn1x_6201`, `e3b_B_D5_un1x_6201` | `phan_tich_doc_lap_u5.py`, `hd_e3b.py` → `analysis/tang_dg/u5v2/n1/phan_tich_n1.py` (no repair) and `sb_chuan_bi/phan_tich_n1_sb.py` (with repair) → `bang_n1.md`, `bang_n1_sb.md`; `external_baselines/cellpack/analysis/phan_tich_cp_n1.py` |
| Sec. 6.6, Table tab:real, Fig. fig:real (five Allen cells; seeds 6541–6543) | `botri_chung_v122u5n14.py` with `specs/n14/dac_ta/*.json` | `u5_n14d5_d5_707681_U5b_6543`, `..._td_6543`, `..._nn1x_6543`, `..._unc8_6543` | `analysis/n14d5/phan_tich_n14d5.py` → `N14D5.json`; render `figures/render/ve_r3.sh` |
| Sec. 6.7, runtime (exclusive re-run, K08 seeds 6511–6513, TV-real) | same commands as above | — | `analysis/do_gio/phan_tich_do_gio_v2.py` (+ `giam_doc_quyen.py` exclusivity monitor, `so_bit.py` bitwise re-run comparison) → `bang_k08_v2.csv`, `bang_tv_v2.csv`, `BANG_GIO_DOC_QUYEN_v2.md` |
| Appendix app:coverage (rule coverage) | — | — | `results/so_lieu/nguon/DoPhu_Luat_BoDung_20261123.md` (+ its review `RaSoat_DoPhu_Luat_20261123.md`) |
| every number printed in the paper | — | — | `results/so_lieu/tao_so_lieu.py` → `so_lieu.tex` |

`results/commands/paper_runs.tsv` lists all 513 runs behind Tables tab:w1, tab:abl, tab:plant, tab:ext, tab:real and the CX comparison (run id, paper item, command). The commands are copied from the run logs; only paths were rewritten to repository-relative paths (`builder/`, `specs/`, `--out runs/<id>`) and the interpreter path to `python`. The plant runs also set `IGL_PHU=<core list>`, which moves `igl.signed_distance` into a helper process (`builder/igl_tho.py`) pinned to fixed CPU cores; the helper evaluates the same function on the same float64/int64 inputs, and the variable is omitted from the listed commands. Development runs (seeds 6501–6505, 6597–6598) are not listed.

## Running one cell

Requirements: Linux, one CUDA GPU with 24 GB (RTX 4090 class), and the environment in `environment/` (`requirements_b2demo.txt`). Run every command from the repository root; the builder writes its caches under `bo_dung_chung/` and `cache/` relative to the working directory.

Place one K08 cell with PA-PLACE (seed 6511; about 11–12 minutes on one RTX 4090 in the exclusive timing run of Sec. 6.7):

```bash
GO_LOI=1 PYTHONNOUSERSITE=1 PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
python builder/botri_chung_v122u5e.py --dac_ta specs/v12/dong_vat_v12_khe0.8.json --seed 6511 --N1 128 --K 4 \
  --thu_vien T --so ti_the=600 --so lysosome=200 --cach tu_dong --tang2 --khong_tu_do --eval_dung_tu_dau --hd_truoc_sau \
  --out runs/u5_r_U5b_k08_6511 --u5 --u5_lambda 1.0 --u5_kappa 0 --u5_chon gc --u5_er --u5_tam
```

PLACE with K = 8 on the same cell: drop the `--u5…` flags and use `--K 8` (run id `u5_r_P8_k08_6511`). The output directory contains `thong_tin.json` (placement report, PROJECT feasibility, per-part pre/post Hausdorff distances in `ban_u5b.hd_truoc_sau`), `mesh_vo_bao.npz` (final meshes) and `mesh_vo_bao_truoc.npz` (pre-projection meshes).

Score the cell with the independent checker:

```bash
python checker/cham_mot.py runs/u5_r_U5b_k08_6511 --vat dong_vat     # plant cells: --vat thuc_vat
```

The verdicts are written to `checker/cham_ra/<name>/kiem3.csv` (one row per rule: pass, fail, partial, not checkable; plus `G_c` and the template deviations `D_i`). The paper scores every method through the same frozen path: a blind copy (`checker/tao_ban_mu.py`), mesh sanitation (`checker/cham_lo_vs.sh` → `ve_sinh_luoi.py`), optionally the uniform repair (`checker/cham_lo_vsb.sh` → `sua_mesh_dn.py`), then `checker/cham_lo.sh` → `cham_mot.py`. The batch scripts check the md5 of the 16 frozen files listed in `DONG_BANG_BO_KIEM_V3.txt` (section A) and expect the original directory layout (a directory named `bo_kiem_v3` next to the GPU scheduler `chay_gpu.sh` of the compute server, which is not part of this repository); `cham_mot.py` runs on its own.

## Regenerating the tables from `results/`

The reviewed result files are in `results/so_lieu/nguon/`; `results/so_lieu/NGUON.md` (Vietnamese) states for each file where it came from, which table it feeds and how it was cross-checked. The paper uses only LaTeX macros generated from these files:

```bash
python3 results/so_lieu/tao_so_lieu.py      # writes results/so_lieu/so_lieu.tex
python3 results/hinh/fig4_paired_seeds.py   # writes the paired-seed figure next to the script
```

Both scripts need only the Python standard library (the figure script also needs matplotlib). The header of `so_lieu.tex` records the md5 of every source file, so a regenerated file can be compared with the shipped one. Recomputing the result files themselves from raw runs requires the per-run outputs (meshes and reports, about 29 GB), which are not in this repository; the analysis scripts in `analysis/` read them from the original directory layout (see below).

## Real cell containers (Allen hiPSC)

The five real containers are hiPSC cells from the Allen Institute for Cell Science hiPSC single-cell image dataset (Viana et al., Nature 2023), used under the Allen Institute terms of use and citation policy. This repository does **not** contain any Allen image, segmentation or mesh. To rebuild the containers:

1. Download, from the Allen Cell Explorer (allencell.org, hiPSC single-cell image dataset), the single-cell segmentation crops (`crop_seg`, nucleus and membrane channels) of the cells with CellId **798836, 989830, 593409, 707681, 997092**. `specs/n14/te_bao/<CellId>/meta.json` records for each cell its structure name, PlateId, FOVId and the `crop_seg` file name.
2. Run `specs/n14/chuan_bi_n14.py`. It builds the membrane and nucleus meshes, the exact signed distance fields on the projection grid (`sdf.npz`, md5 recorded as `md5_sdf` in `meta.json`) and the nucleus template (`khuon_nhan.npz`), and writes the specifications. The selection rule (stratified by cell size, seed 20261119) and the eligibility table of the 200 candidate cells are in the script header and in `specs/n14/chon_te_bao.csv`; `specs/n14/chon.json` records the five chosen cells. The script reads the crops and the 200-cell sample table from the original server layout (variable `ALLEN`) and the cached membrane field of the synthetic animal cell written by the builder (`SDF_D4`).
3. The specifications `specs/n14/dac_ta/dong_vat_n14_<CellId>.json` point to `sdf.npz` and `khuon_nhan.npz` with absolute paths of the compute server (fields `vo.tep_sdf`, `khuon[…].tep`, `bo_phan[0].co_dinh_tep_sdf`). They are shipped unchanged to keep their md5; edit these three fields to your local files before running.
4. Organelle counts at D5-like density: `specs/n14/d5_density/so_luong_n14d5.json` (from `kha_thi_n14d5.py`).

## External baselines

* **cellPACK** (`external_baselines/cellpack/`): `bridge/cau_noi_cp_e3.py` exports each cell to a cellPACK recipe (`xuat`), `runner/chay_cp_gpu.py` runs cellPACK, `bridge/cau_noi_cp_e3.py chon` / `nhap` select and import the poses, and `botri_chung_v122e3.py` then applies the same PROJECT, repair and checker as for every method. The queue script `bridge/hang_n1.sh` records the exact sequence and arguments (cellPACK seed 10·seed + k, `--thu_tu ban_kinh_giam --bo tac_gia`). The GPU port of cellPACK is only a patch: `gpu_patch/v4.diff` is the `git diff` against `mesoscope/cellpack` commit `a576c5c1e81a16d7b908bfba9168e708f4730ce7`, and `gpu_patch/v4_gpu_tang_toc.py` is the one new module (`cellpack/autopack/gpu_tang_toc.py`); `gpu_patch/v4.md5` lists the md5 of the patched files. In the queue, poses of the patched version are imported only after an equivalence gate passes: selected cells (chosen by `bridge/chon_cap_e3.py`) are run with both the original and the patched cellPACK and compared bit for bit by `runner/so_sanh.py` (placed parts, positions, rotations, grid fingerprints, random-generator state), and every pair must be identical. cellPACK is distributed under the MIT License, Copyright (c) 2021 Graham Johnson, Ludovic Autin, Megan Riel-Mehan; the full notice is in `external_baselines/cellpack/LICENSE.cellpack` and applies to the cellPACK code contained in the patch.
* **Packmol** (`external_baselines/packmol/`): `cau_noi_pk_e3.py` (export, select, import), `pk_dau_vao.py` (writes the Packmol input), `pk_doc.py` (reads the Packmol output); `hang_pk_n1.sh` records the sequence (Packmol 100.0.0, `packmol -t 1`, default author settings, Packmol seed 10·seed + k). Packmol itself (MIT License, L. Martínez and contributors) is not redistributed; `SHA256SUMS_packmol-100.0.0.txt` lists the checksums of the official release files.

## Limitations of this release

* A CUDA GPU is required: the builder runs PLACE and PROJECT in PyTorch on the GPU (24 GB was enough for every run in the paper). The checker and the analysis scripts run on the CPU.
* Byte-identical copies keep the absolute paths of the compute server in several files (`~/Vision_team_data_189/...`, `/home/<user>/...`). They matter for the N14 specifications (see above), for the checker's AI-template lookup (`checker/do_khuon.py`, variable `GOC_SV`), for the batch scoring scripts, and for the analysis scripts, which read per-run outputs from the original layout. `SOURCE_PATHS.tsv` gives the original location of every file, so the layout can be recreated.
* The library of AI-generated templates for D4A2 is included in `data/ai/khuon_chon/` (see the section on shipped data below). The checker's plant templates (`ai/khuon_tv/`) are mesh files and are not included. Procedural templates (`--thu_vien T`) need nothing extra.
* The cellPACK and Packmol bridges check the md5 of the cached membrane field of the synthetic animal cell (`bo_dung_chung/v121/sdf_mang/sdf_75005fbc7d88018f.npz`) before they run, and `specs/n14/chuan_bi_n14.py` reads the same file; the builder writes it on its first animal run, and a different GPU or library version may produce a field with a different md5.
* Per-run outputs (meshes, reports, checker outputs) are not included because of their size; `results/` holds the reviewed per-seed and aggregate numbers.
* The CX baseline is not bitwise deterministic across re-runs (non-deterministic GPU atomic sums feed back into its periodic approximate projection; Sec. 6.7); PA-PLACE, P4 and P8 reproduced placements and meshes exactly.

## License and citation

The code in this repository is released under the MIT License (`LICENSE`). The cellPACK patch also contains cellPACK code under its own MIT License (`external_baselines/cellpack/LICENSE.cellpack`). Please cite the paper if you use this code (`CITATION.cff`); the DOI will be added after publication.

## Data shipped with the code, and paths of the original machine

- `data/ai/khuon_chon/` is the library of AI-generated templates used by the D4A2 configuration (`--thu_vien A`), with its index files `danh_sach.json` and `danh_sach_v2.json`. The builders expect it at `ai/khuon_chon/` relative to the working directory; copy or link `data/ai` there before running.
- The membrane SDF cache `bo_dung_chung/v121/sdf_mang/sdf_75005fbc7d88018f.npz` (112 MB, md5 `ee14b688e6225311735c3ebfbeecc902`) is not included because of its size. The builder writes it on the first run from the cell specification; the cellPACK and Packmol bridges check its md5 before they start.
- The data of the real Allen hiPSC cells are not redistributed; see `specs/n14/` for the cell ids and the preparation script, and obtain the data from the Allen Cell Explorer under its terms of use.
- The final meshes shown in the figures are available from the corresponding author on request.
- The scripts are the versions that produced the paper's numbers (comments removed), so they still contain absolute paths of the original machine (for example `/home/<user>/...` and `~/Vision_team_data_189/...`). Adjust these paths, or create matching links, before running them elsewhere; the Python environment is listed in `environment/`.
