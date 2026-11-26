# HINH_CAD — hình 3D cho bài CAD (phiên thực nghiệm, 24/11/2026)

Thư mục server: `~/Vision_team_data_189/sonhh20/ch3_tebao_dep/bo_dung_chung/hinh_cad/`. Bản sao trên Mac: `Inside3D/BaiBao_Bai1_CAD/hinh/`.
Mọi render: Blender 4.2.1 LTS, Cycles GPU, chạy qua `CHUNG=1 chay_gpu.sh auto <GB>` và `taskset -c 30,31`. Đỉnh GPU đo theo PID (`scripts/do_gpu_mem.sh`): R1 2 936 MiB (đặt chỗ 4 GB), R3 2 826 MiB (4 GB), R2 4 522 MiB (5 GB).
Ảnh Blender có nền trong suốt; `scripts/ghep.py` ghép lên nền trắng.
Trạng thái số liệu trong hình: **CHƯA rà soát chéo** (HD_i đã khớp bộ dựng tới 5e-5 µm, xem dưới).

## R1 — `R1_k08_P8_vs_PAPLACE.png` (3310 × 2648 px)

**Nguồn (cấu hình k08, seed 6512):**

| nhánh | npz sau chiếu (md5) | npz trước chiếu (md5) |
|---|---|---|
| P8 | `tang_dg/kq/u5_r_P8_k08_6512/mesh_vo_bao.npz` (44e6312714aa4af7fe8196ceb92723d5) | `mesh_vo_bao_truoc.npz` (b9c99b9d394c066377b968ece2821b35) |
| PA-PLACE (U5b) | `tang_dg/kq/u5_r_U5b_k08_6512/mesh_vo_bao.npz` (3b396f6e683113f1ada0f4ca77033a78) | `mesh_vo_bao_truoc.npz` (6f844e9c7a9133dc7ecb82af4e4e0de9) |

**HD_i theo bộ phận.** Thư mục `u5_r_*` không có `n2.json`, nên HD_i được tính lại bằng `scripts/hd_n2_dv.py` (md5 5eb72ce9ddc8bed0bb4a80c34de503ab).
Script này exec đúng hai hàm `mesh` và `kc` của thước N2 `tang_dg/u5tv/n2_tv.py` (md5 14b46b125905736cb4217fd7d7b27a51, có kiểm md5 trước khi chạy).
- Định nghĩa: HD_i = max(max kc(đỉnh sau → mặt trước), max kc(đỉnh trước → mặt sau)), với kc = `sqrt(igl.point_mesh_squared_distance)`.
- Hai mesh trùng hệt thì HD_i = 0.
- Kiểm chéo với `thong_tin.json` → `ban_u5b.hd_truoc_sau` của bộ dựng: lệch tối đa 5e-5 µm ở cả 20 việc k08 (do bộ dựng làm tròn 4 chữ số).
- Kết quả nằm ở `du_lieu/hd_u5_r_{P8,U5b}_k08_<s>.json` (md5 seed 6512: P8 9f9c4fa6…, U5b 9de815fd…).
- "HD max tế bào" = max HD_i trên mọi bộ phận, trừ màng (HD = 0) và tầng 2 (trung_tu, tui_golgi: bộ dựng chép mesh sau vào bản trước, nên không có bản trước thật). Định nghĩa này trùng `hd_ts_tat_ca.max` của `u5v2/doc_lap`.

**Chọn seed: 6512 (không chọn 6511).** Hiệu HD max (P8 − U5b) theo seed, tính từ hd_n2_dv:

| seed | P8 | U5b | hiệu | cách trung vị |
|---|---|---|---|---|
| 6512 | 0.2771 | 0.1838 | 0.0933 | 0.0068 |
| 6518 | 0.2741 | 0.1671 | 0.1069 | 0.0068 |
| 6515 | 0.3595 | 0.2464 | 0.1131 | 0.0129 |
| 6517 | 0.2841 | 0.2079 | 0.0762 | 0.0239 |
| 6519 | 0.3416 | 0.2090 | 0.1326 | 0.0325 |
| 6513 | 0.2805 | 0.2138 | 0.0667 | 0.0334 |
| 6520 | 0.2686 | 0.2048 | 0.0638 | 0.0363 |
| 6516 | 0.3062 | 0.2429 | 0.0633 | 0.0369 |
| 6514 | 0.3556 | 0.1958 | 0.1598 | 0.0597 |
| 6511 | 0.4034 | 0.1864 | 0.2171 | 0.1169 |

- Trung vị của hiệu là 0.1001 µm.
- Seed 6511 có hiệu LỚN NHẤT (cặp có lợi nhất cho PA-PLACE), nên không gần trung vị và bị loại.
- 6512 và 6518 cách trung vị bằng nhau (số chẵn seed, trung vị nằm giữa hai seed này). Chọn 6512 theo hai luật phá hoà, cả hai đều cho 6512:
  - (a) HD max từng nhánh gần trung vị của nhánh hơn: trung vị P8 0.2951, U5b 0.2063. Tổng độ lệch của 6512 là 0.0405, của 6518 là 0.0603.
  - (b) Hiệu nhỏ hơn, tức thận trọng, ít có lợi cho PA-PLACE.

**Bố cục.**
- Cột trái P8, cột phải PA-PLACE.
- Hàng 1: mặt cắt nửa.
- Hàng 2: ba ô phóng.
- Dưới cùng: thanh màu chung.

**Camera (cả hai cột, mọi ô giống nhau).**
- Trực giao, nhìn theo +y, nâng ELEV = 15°, AZIM = 0, khoảng cách 60 µm.
- Mặt cắt: tâm (0, 0, 0), ortho_scale 22 µm, 1600 × 1600 px, 256 mẫu, khử nhiễu OptiX.
- Ghép: cắt theo hộp bao CHUNG của hai ảnh, nên tỉ lệ hai cột giống nhau.

**Mặt cắt nửa** (như `n14/ve_n14_blender.py`):
- Mesh lớn (màng, nhân, ER, golgi) cắt trong shader ở y < 0.
- Bào quan nhỏ giữ nguyên vẹn nếu trọng tâm đỉnh có y ≥ 0.
- Màng alpha 0.10. ER alpha 0.35. Các bộ phận khác đặc.
- Hình học không làm trơn, chỉ smooth shading.

**Ô phóng.**
- Hộp lập phương nửa cạnh 1.2 µm, ortho_scale 2.76 µm, 900 px.
- Tâm hộp là trọng tâm đỉnh mesh sau của ba bào quan nhỏ (ti_the/lysosome) có HD_i lớn nhất của P8, trong số bào quan được giữ ở mặt cắt nửa (y ≥ 0):
  - (1) ti_the_390: HD_i 0.277, tâm (0.57, 5.83, 6.67);
  - (2) ti_the_295: HD_i 0.249, tâm (−4.84, 5.72, 4.59);
  - (3) ti_the_11: HD_i 0.248, tâm (5.56, 2.75, −6.47).
- Hạng 2 toàn tế bào là ti_the_588 (0.260). Nó nằm ở nửa trước (y = −0.84), đã bị bỏ khỏi mặt cắt nên không đánh dấu được, vì vậy không chọn.
- PA-PLACE: cùng ba hộp (cùng vùng không gian). Tên bộ phận giữa hai nhánh KHÔNG tương ứng: trọng tâm cùng tên cách nhau trung vị 11 µm.
- Mỗi ô vẽ mọi ti_the/lysosome có trọng tâm trong hộp: mesh sau đặc, tô theo HD_i. Mesh trước chiếu của chính chúng vẽ bằng khung dây xanh #00A5D9, dày 0.008 µm, không đổ bóng. Phần dây nằm ngoài khối đặc là phần bị gọt.
- Không vẽ mặt mờ của mesh trước: ở vùng không bị gọt nó trùng mặt với mesh sau, gây nhiễu z.
- Mesh lớn (nhân, ER, màng) không vẽ trong ô phóng.
- Nhãn ô ghi max HD_i trong hộp: P8 0.277 / 0.249 / 0.248; PA-PLACE 0.143 / 0.114 / 0.080 µm.
- Vòng tròn đánh số ở hàng 1 là vị trí tâm hộp, chiếu trực giao.

**Thang màu.**
- Tuần tự `magma_r` của matplotlib, lấy đoạn [0.04, 0.80] (bỏ đầu gần trắng và đuôi gần đen), 256 mức, trong `du_lieu/lut_magma_r.json` (md5 476fcbfd…).
- Tuyến tính từ 0 tới VMAX = 0.2771 µm, là max HD_i của cả hai tế bào (đạt ở ti_the_390 của P8). Một thang chung cho cả hai cột.
- Bộ phận tầng 2 và nhân con có HD_i = 0 nên tô màu đầu thang. Riêng tầng 2, giá trị 0 nghĩa là "không có bản trước thật", không phải "không méo".

**Tái lập.**
1. `bash scripts/ve_r1.sh` (md5 d7b24688…), chạy qua `CHUNG=1 chay_gpu.sh auto 4`.
2. `python scripts/cfg_r1.py` rồi `python scripts/ghep.py du_lieu/cfg_r1.json`.
- Script Blender: `scripts/ve_hd_blender.py` (bản cuối md5 e0006397f872529dc0e0bc04ac394d1a, dùng chung cho R1–R3; R1 đã render lại bằng bản này).
- File `.blend` của hai mặt cắt: `r1/cat_P8.blend`, `r1/cat_U5b.blend`.
- Ảnh thành phần: `r1/*.png`. Ảnh ghép md5 00ab291dbb261f6eec7941791cb2bb89.

**Ghi chú khi đọc hình.**
- Ở P8, nhân (HD_i 0.188) và golgi (0.211) cũng bị méo rõ hơn PA-PLACE (0.101 và 0.184). Đây là số liệu thật, không phải lỗi màu.
- Mesh trích ở bước lưới h = 0.18 µm nên ô phóng thấy rõ mặt tam giác.

## R3 — `R3_N14_707681_PAPLACE_vs_PLACE.png` (3290 × 2689 px)

**Chọn tế bào và seed.**
- Biến thể d5 (D5-like), nhánh PA-PLACE (U5b). 15 cặp (tế bào, seed).
- Thước là `HD_max_dat` (= `u5v2/doc_lap/<id>.json` → `hd_ts_tu_the.max`, thước chính của bảng N14D5): HD max trên các bộ phận được ĐẶT, không gồm ER (ER bị màng thật cắt) và không gồm nhân cố định.
- Trung vị 15 giá trị = 0.1047 µm, đạt đúng ở **tế bào 707681, seed 6543**. Không phải chọn hình đẹp.
- Cột PLACE (td) dùng cùng tế bào, cùng seed: HD_max_dat = 0.2062 µm.

**Nguồn.**

| nhánh | npz sau chiếu (md5) | thong_tin.json (md5) |
|---|---|---|
| PA-PLACE (U5b) | `tang_dg/kq/u5_n14d5_d5_707681_U5b_6543/mesh_vo_bao.npz` (f96bc0075ba98efdf8d0eebf41c0ec38) | 9f3affff1b3a2efb9b327e460bf497bc |
| PLACE (td) | `tang_dg/kq/u5_n14d5_d5_707681_td_6543/mesh_vo_bao.npz` (87dcda5bc7bffa4c4264f675fdf752c8) | 1dd56bb16f78f34949f3c8c285c04faf |

**HD_i.**
- Việc N14 KHÔNG giữ mesh trước chiếu, nên không chạy lại được hàm N2.
- HD_i lấy từ bộ dựng (`thong_tin.json` → `ban_u5b.hd_truoc_sau`), qua `scripts/hd_tu_bo_dung.py` (md5 6df1817f…).
- Định nghĩa trùng N2. Trên 20 việc k08, giá trị bộ dựng khớp `hd_n2_dv.py` tới 5e-5 µm.

**Bố cục.**
- Cột trái PA-PLACE, cột phải PLACE.
- Hàng 1: nhìn 3D toàn tế bào, không cắt, nâng 35°.
- Hàng 2: mặt cắt nằm ngang, bỏ nửa trên z > 0, nâng 55°. Tế bào thật dẹt (z trong ±6.5 µm), nên cắt theo y sẽ bỏ gần hết bào quan.

**Camera.**
- Trực giao, AZIM 0, nhìn theo +y, ortho_scale 27 µm, 1600 × 1200 px, 256 mẫu.
- Mỗi hàng cắt theo hộp bao chung của hai cột, nên tỉ lệ hai cột giống nhau.

**Màu.**
- Màng thật #D9D9D9, alpha 0.12. Nhân thật #B8B8D0, đặc. ER #8FA3B5, alpha 0.30. Ba bộ phận này KHÔNG tô theo HD.
  - Lý do với ER: HD của ER ≈ 3.6 µm vì bị màng thật cắt, sẽ làm hỏng thang màu.
- Bộ phận được đặt (ti thể, lysosome, golgi, trung thể, nhân con) tô theo HD_i.
- Thang `magma_r` [0.04, 0.80], từ 0 tới VMAX = 0.2062 µm (max HD_i của bộ phận đặt ở cả hai tế bào).
- Tầng 2 (trung_tu, tui_golgi) không có trong hd_truoc_sau, nên tô màu đầu thang.

**Tái lập.**
- `bash scripts/ve_r3.sh` (md5 44781c8b…) qua `CHUNG=1 chay_gpu.sh auto 4`. Đỉnh GPU đo được 2 826 MiB.
- Sau đó `python scripts/ghep_luoi.py du_lieu/cfg_r3.json` (md5 ghep_luoi bec6cdb0…).
- File `.blend`: `r3/{3d,cat}_{U5b,td}.blend`. Ảnh ghép md5 3f32b9dab7ecce1b11f4bc5ff7f285f0.

## R2 — `R2_TVreal_7311_ba_nhanh.png` (3230 × 4639 px)

**Nguồn: TV-real, seed 7311, ba nhánh còn giữ npz.**

| nhánh (nhãn hình) | npz sau chiếu (md5) | n2.json (md5) | HD_i rắn max |
|---|---|---|---|
| PA-PLACE+F (`u5ktt`) | `tang_dg/kq/u5_tv_c_real_u5ktt_7311/mesh_vo_bao.npz` (7b49a5d2c50b237a189151baa8a800b4) | 27bef59e55a34698f0bf7665024d0773 | 0.069 |
| PLACE (`td`) | `tang_dg/kq/u5_tv_c_real_td_7311/mesh_vo_bao.npz` (24f2ae0bf96f385568bfaa7d9fcd8182) | 110c48d925c63eec12a682fb57c79e77 | 0.414 |
| thử–loại (`nn1x`) | `tang_dg/kq/u5_tv_c_real_nn1x_7311/mesh_vo_bao.npz` (eb2c5294dbe2d13e97eba00aa64cb43d) | 44b198726bd0fa623ea0c97c57a60ea7 | 1.248 |

Nhãn hàng thứ ba trong hình là "Rejection". Cần đối chiếu với tên nhánh trong bài: bảng thực vật dùng "Rej.+soft", có thể là `nn1xm`, còn chỗ giữ hình ghi `nn1x`.

**Cột SA (`un1x`): KHÔNG thêm được.**
- `kq/u5_tv_c_real_un1x_7311/` chỉ còn `n2.json`, `thong_tin.json`, `trang_thai_luat.json`.
- Lô chấm `u5tv/cham/lo_10`, mã mù 58cd5ef907d9 (sau) / 014ecc1fc40a (trước): thư mục `vs/vs/…` và `vs/sb/…` chỉ còn json và dấu, không còn `mesh_vo_bao.npz`.
- `bo_kiem_v3/cham_ra/vsu5tv10_mu_sau__58cd5ef907d9/mesh.npz` là liên kết gãy.
- Tìm trên cả ổ đĩa không thấy bản sao nào khác.
- Mesh sau chiếu là MC của trường sau phép chiếu, nên không dựng lại được nếu không chạy lại bộ dựng.

**HD_i.**
- Lấy từ `n2.json` của chính việc (thước N2 `n2_tv.py`, md5 14b46b12…), qua `scripts/hd_tu_n2.py` (md5 c8fe550c…).
- Bộ phận rắn (lục lạp, ti thể, peroxisome, golgi, nhân): `bo_phan[k].hd`.
- ER và không bào lấy từ `vung`, nhưng không dùng để tô màu.

**Bố cục.**
- Ba hàng (một nhánh mỗi hàng), cùng camera.
- Mặt cắt nửa: bỏ y < 0. Vách, màng, không bào cắt trong shader. Bào quan nhỏ giữ nếu trọng tâm y ≥ 0.
- Vách #C8C8C8, alpha 0.08. Màng #B0C4DE, alpha 0.08. Không bào #9FC5E8, alpha 0.15. Ba bộ phận này không tô HD.
- ER KHÔNG vẽ: lưới tấm mỏng 3.6 triệu mặt che kín lớp tế bào chất ở bản thử.
- Bộ phận rắn đặc, tô HD_i.

**Camera.**
- Trực giao, nhìn theo +y, nâng 15°, ortho_scale 84 µm (tế bào 80 × 32 × 32 µm).
- 3000 × 1500 px, 256 mẫu, khử nhiễu OptiX.

**Thang màu.**
- `magma_r` [0.04, 0.80], từ 0 tới VMAX = 1.2479 µm (max HD_i rắn của ba nhánh, đạt ở nhánh nn1x). Một thang chung.
- Hệ quả: nhánh PA-PLACE+F (max 0.069) gần như toàn màu đầu thang. Đó là số liệu, không phải lỗi màu.

**Tái lập.**
- `bash scripts/ve_r2.sh` (md5 e3d885ab…) qua `CHUNG=1 chay_gpu.sh auto 5`.
- Sau đó `python scripts/ghep_luoi.py du_lieu/cfg_r2.json`.
- Không lưu `.blend`: mỗi cảnh khoảng 31 triệu tam giác, `.blend` khoảng 1 GB. Bật `BLEND_ON=1` nếu cần.
- Ảnh ghép md5 4f1b421c8a6b365397dfa3cf5a7a6332.

## Bản in v2 (24/11 tối): chỉ ghép lại, KHÔNG render lại

Theo yêu cầu phiên viết:
- Chữ ≥ 7 pt khi in rộng 16 cm. Ảnh rộng 3150 px ở 500 dpi, nên 7 pt ≈ 49 px. Cỡ chữ dùng ở đây:
  - tiêu đề 84 px;
  - phụ đề, nhãn hàng, nhãn thang 60–62 px;
  - vạch thang, nhãn ô phóng 56 px;
  - ghi chú 54 px;
  - nhãn ô (a)… 66 px.
- Chính tả Mỹ.
- Thêm nhãn ô (a)(b)…
- R2: bỏ tiêu đề và ghi chú SA. Nhãn hàng là "PA-PLACE+F", "PLACE", "Rejection".
  - Ghi chú: PLACE (td) và PA-PLACE chạy với --K 4. Nhãn không ghi K.
- Xuất PNG và TIFF LZW, 500 dpi, khổ 16.0 cm.

Script `scripts/ghep_v2.py` (md5 b35399b4…), cấu hình `du_lieu/cfg_r{1,2,3}_v2.json`. Các ảnh thành phần vẫn là `r1/`, `r2/`, `r3/` như trên.

| hình | PNG (md5) | TIFF (md5) | kích thước |
|---|---|---|---|
| R1 | `R1_K08_P8_vs_PAPLACE_v2.png` (e5ce8c1373b7bbb882274ce419e4b27d) | a5abf2959eff110222b9ae568b264dbe | 3150 × 2686 (16.0 × 13.6 cm) |
| R2 | `R2_TVreal_7311_v2.png` (9201bb705f95cc0ce6639bde250f573a) | 3d7473a7070addb136b3af13c1961088 | 3150 × 4338 (16.0 × 22.0 cm) |
| R3 | `R3_N14_707681_v2.png` (0293babb1c4f2c134b7b0a84857e50c0) | 6e0518c6782612ee4eddf6b2f87b18b0 | 3150 × 2663 (16.0 × 13.5 cm) |

Chú thích bài cần tự nêu, vì đã bỏ khỏi hình:
- R2: ER không vẽ; thang màu chung 0–1.248 µm; không có nhánh SA (mesh seed 7311 đã bị xoá).
- R1: ô phóng ghi max HD_i trong hộp.

## Bản in v3 (24/11 tối): ghép lại, không render lại

Lý do: bài đặt hình ở \linewidth ≈ 13.7 cm, nên với ảnh rộng 3150 px thì 1 pt ≈ 8.1 px.

Thay đổi so với v2:
- Mọi cỡ chữ ×1.2 (khoá `"co_chu": 1.2` trong cấu hình). Cỡ nhỏ nhất 65 px ≈ 8.0 pt khi in 13.7 cm; vạch thang và nhãn ô phóng 67 px ≈ 8.3 pt.
- R1 bỏ dòng ghi chú cuối. Nhãn ô phóng chuyển xuống góc dưới phải (`"nhan_o_duoi": true`) để không đè nhãn (c)/(d).

Script `scripts/ghep_v2.py` (md5 b0e5b6a56d320f4515a7ed7307c086fc). Khi không có `co_chu` và `nhan_o_duoi`, script cho ra đúng bản v2. Cấu hình `du_lieu/cfg_r{1,2,3}_v3.json`.

| hình | PNG (md5) | TIFF LZW (md5) | px |
|---|---|---|---|
| R1 | `R1_K08_P8_vs_PAPLACE_v3.png` (8a4118d5ccb39bcd38a2631b132bf254) | 895a113ec12b1c6a53694855e6a88d78 | 3150 × 2578 |
| R2 | `R2_TVreal_7311_v3.png` (fc76a6acbb5c5db1af5b080097329ff3) | 01c23c09f7d59ebcde22f644c8285d03 | 3150 × 4334 |
| R3 | `R3_N14_707681_v3.png` (575763d1ea20f0e04fc186af4fd289fb) | b88ac716c27206c7db49e9b44fb08395 | 3150 × 2706 |

Chú thích R1 cần nêu thay cho dòng ghi chú đã bỏ: K08, seed 6512; nhãn ô phóng là max HD_i trong hộp (µm); khung dây xanh là mesh trước chiếu.
