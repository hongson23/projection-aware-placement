# Nguồn số liệu của bài CAD (sinh macro bằng `tao_so_lieu.py` → `so_lieu.tex`)

Quy tắc: `main.tex` và `main_vi.tex` chỉ dùng macro. Muốn đổi số thì đổi tệp nguồn rồi chạy lại `python3 so_lieu/tao_so_lieu.py`.
md5 của mọi tệp nguồn được ghi ở đầu `so_lieu.tex` mỗi lần sinh.

| tệp (nguon/) | lấy từ server (`…/ch3_tebao_dep/bo_dung_chung/tang_dg/`) | dùng cho | trạng thái rà soát |
|---|---|---|---|
| `ket_qua_rong.json`, `ket_qua_rong.md` | `u5v2/` | W1: U5b, P4, P8, U5b8, a≡1 (λ=1), λ=0 ở k08 / D5 / D4A2, seed 6511–6520 | ĐÃ rà độc lập: số tính lại khớp `RS_U5R_du10seed_20261121.md` (21/11 22:49) và `Inside3D/RaSoat_TongHop_KetQua_20261122.md` |
| `RS_U5R_du10seed_20261121.md` | Mac `ThuNghiem/ra_soat_u5r/ket_qua/` | đối chiếu | bản tính lại của người rà |
| `ket_qua_u5c.json`, `ket_qua_u5c.md` | `u5v2/u5c/` | đối chứng cơ chế k08: a≡1 (λ 0,3/1/3), a≡½, W∝V | bảng 10 seed đã tính lại độc lập 23/11 (`BaiBao_Bai1_CAD/RA_SOAT_BAN1.md` §0, S6) |
| `ket_qua_cx.json`, `ket_qua_cx.md`, `chon_dev.json` | `u5x/` | CX vòng 1 (CXe = BL1l1m800) | ĐÃ rà chéo: `ThuNghiem/doi_chung_cx_20261121/RA_SOAT_CX_CHINH_20261122.md` (khớp từng ô) |
| `chon_dev2.json`, `chon_dev2.md` | `u5x/` | CX vòng 2 (chọn cấu hình dev) | CHƯA rà; lượt chính CX2 → `ket_qua_cx2.json` (TBD) |
| `so_lieu_vi.tex` | sinh cùng `so_lieu.tex` | bản Việt: dấu phẩy thập phân | — |
| `tay.json` | Mac, chép từ tệp rà soát | phần vi phạm còn lại k08, tổng bộ phận <50 %, D5 chậm 11 %, CX best bất khả 6–846, tham số thiết kế (h, số bào quan, khe) | mỗi mục ghi tệp + mục nguồn; chỉ số đã rà |

**Chưa có (macro `\TBD`):** CX vòng 2 lượt chính; thực vật lượt chính (u5tv); bảng thống nhất N1 (cellPACK, Packmol); N14 mật độ cao; đo giờ độc quyền rút gọn.

**Định nghĩa trong macro:**
- `…Rel` = |trung vị hiệu tương đối theo từng seed| (%), không phải tỉ số trung vị.
- `…Win` = thắng/thua theo cặp seed (hoà ghi trong chú thích của macro).
- `HD p95` ở bảng W1 = `HD_tu_p95` (trừ ER); ở đoạn CX = `HD_p95_tat_ca` (mọi bộ phận). Bài đã ghi rõ hai định nghĩa.

| `ket_qua_cx2.json`, `ket_qua_cx2.md` (md5 5714c786) | `u5x/` | CX ba bản (CXe, CX2e, CX2b) so với U5b, P4, P8 — Bảng tab:cx | ĐÃ rà chéo: `ThuNghiem/doi_chung_cx_20261121/RA_SOAT_CX2_20261123.md` (khớp từng ô; câu được phép §6; giới hạn §5 a–f). Số Σ d50, ngân sách, thời gian đặt, số lần chiếu lấy từ bản rà đó (tay.json, tiền tố y…). Rút lại "6–846": đúng là 2–846. |
| `DoPhu_Luat_BoDung_20261123.md`, `RaSoat_DoPhu_Luat_20261123.md` | Mac `BaiBao_Bai1_CAD/` (chép vào nguon/) | Bảng độ phủ luật (phụ lục `app:coverage`), tiền tố macro `cv…` trong tay.json | Bảng DoPhu CHƯA rà riêng; bản rà RaSoat_DoPhu (một vòng, đếm lại từ dữ liệu gốc) — chỉ dùng số đã được bản rà xác nhận/sửa (4.3b 0,107–0,191; M.c₀ 100 % D4A2, không gọi lỗi bộ dựng; N14 45 thật + 45 tổng hợp; 2.7c không thuần do lưới) |

**Hình (24/11):** `hinh/cap_seed_k08.pdf` (Hình cặp seed K08) do `hinh/ve_cap_seed_k08.py` vẽ, đọc DUY NHẤT `nguon/ket_qua_cx2.json` (md5 5714c786, ĐÃ rà chéo: RA_SOAT_CX2_20261123.md). Giá trị theo seed: PA-PLACE = phía A của cặp `U5b-CX*`; P4, P8 = phía B của `CXe-P4`, `CXe-P8`; trung vị khớp macro (HD max 0.20635 / 0.29515).
**Bảng thực vật (tab:plant):** macro `\tv<Mức><Nhánh><Thước>` đang là TBD (sinh trong `tao_so_lieu.py`); khi có `nguon/ket_qua_tv.json` ĐÃ rà thì điền.

**Hình render 3D (24/11, phiên thực nghiệm; mô tả đầy đủ: `hinh/HINH_CAD.md`).** Bản PNG gốc → JPEG q93 dùng trong bài:
| hình | tệp bài | PNG gốc (md5) | số in trong hình | trạng thái rà |
|---|---|---|---|---|
| Hình 1 (fig:k08) | `hinh/R1_k08.jpg` | `R1_k08_P8_vs_PAPLACE.png` (00ab291d) | HD max tế bào P8 0,277 / PA-PLACE 0,184 (seed 6512) | khớp giá trị theo seed của `ket_qua_cx2.json` (ĐÃ rà): P8 0,2771, U5b 0,1838. HD_i từng bộ phận khớp bộ dựng tới 5e-5 µm; bản thân hình CHƯA rà |
| fig:plant | `hinh/R2_tv_real.jpg` | `R2_TVreal_7311_ba_nhanh.png` (4f1b421c) | max HD_i 0,069 / 0,414 / 1,248 | CHƯA rà (lượt thực vật chưa rà); chú thích có \TBD |
| fig:real | `hinh/R3_n14.jpg` | `R3_N14_707681_PAPLACE_vs_PLACE.png` (3f32b9da) | max HD_i 0,105 / 0,206 | CHƯA rà (N14 chưa rà); chú thích có \TBD |

| `N14D5.json` (md5 6cd2fa39; từ `n14d5/ket_qua/`) | `../n14d5/ket_qua/` | Bảng tab:real, macro `\nf…` | ĐÃ rà chéo (`ThuNghiem/RaSoat_KetQua_N14D5_20261124.md`, TongKet §21): mọi số khớp. Câu được phép: "U5b ≥ PLACE ≫ SA, thử–loại". HD max U5b so với PLACE chỉ "không kém" (thước độc lập 3/2/0 tế bào) ⇒ KHÔNG in thắng/thua HD max với PLACE. Ngân sách thử–loại 0,21×, SA 0,52× (tay.json `nfBudRej`, `nfBudSA`); ER 45 % ngoài màng từ đầu, mất ~81 % mọi nhánh |

| `ket_qua_tv.json` / `.md` (md5 aedc38f5; = u5tv/KET_QUA_DAY_DU, phan_tich_tv_u5.py 304bf2f0) | `u5tv/` | Bảng tab:plant, macro `\tv…` (tính trung vị và thắng/thua theo seed từ JSON từng việc; đã đối chiếu với bảng .md) | ĐÃ rà chéo (TongKet §25; `ThuNghiem/u5tv_20261122/RA_SOAT_MOC_1_3_20261124.md` mục "Đủ 6 seed"). Câu được phép: TV-D5 PA-PLACE < PLACE, SA ở 6/6 seed, < thử–loại ở 5/6; TV-real < cả ba ở 6/6; TV-D1 chỉ "no harm"; KTT nhận D1 3/3, D5 0/6, real 0/6 (không viết "nhờ KTT"); câu đánh đổi 2.7b/2.7f đã rà. Bài không dùng ngưỡng "rõ" (chỉ nêu số seed), nên không cần khai ngưỡng 6/6 định sau mốc 2/3. |
| `bang_n1.md` (ba8bcc60; 25/11 01:32), `bang_n1_sb.md` (50192b2b; 02:20) | `u5v2/n1/`, `u5v2/n1/sb_chuan_bi/` | Bảng tab:ext, macro `\nl…` (không sửa mesh), `\nlsb…` (có sửa mesh), tay.json `nl…` | ĐÃ rà chéo (`ThuNghiem/RaSoat_BangN1_20261125.md` + `_tinhlai.json`; TongKet §26, §28). Đã khai: cellPACK bỏ khoảng 27 % bộ phận ở D5; Packmol không hội tụ ở D5; PROJECT gốc cho nhánh ngoài; PLACE D5 HD n = 9; 7 cấu hình xuất hỏng và 1 lượt cellPACK hết RAM; M.c0 đổi ở 35 lượt D4A2 khi có sửa mesh (ưu thế luật của U5b ở D4A2 là hiện vật lưới). |

| `bang_k08_v2.csv`, `bang_tv_v2.csv`, `BANG_GIO_DOC_QUYEN_v2.md` (phan_tich_do_gio_v2.py) | `bo_dung_chung/do_gio/` | Mục 6.7, macro `\tm…` (trung vị pp theo nhánh; trung vị tỉ số theo seed so với P8) | ĐÃ rà chéo (`ThuNghiem/do_gio_20261122/RA_SOAT_KET_QUA_DG_20261125.md`, TongKet §29). Câu EN/VI dùng nguyên văn §5. Đã bỏ câu "changes below Y" (chưa chấm lại 2 mesh CX phân kỳ). |
