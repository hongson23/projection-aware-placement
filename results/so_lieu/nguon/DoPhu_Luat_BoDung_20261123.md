# Độ phủ luật × bộ dựng cuối (Bài 1 ch3, bài CAD) — 23/11/2026

Trạng thái: CHƯA rà soát chéo (bảng dựng từ mã + đặc tả, không dựa tài liệu nội bộ).

Nguồn (server, chỉ đọc): `~/Vision_team_data_189/sonhh20/ch3_tebao_dep/`
- Bộ kiểm: `bo_kiem_v3/kiem3.py` (đăng ký luật bằng `B.them(...)`), ngưỡng sự thật tách từ `so_trong_van_ban` qua `bo_kiem_v3/luat_v12.py`, tham số `bo_kiem_v3/luat_v12/tham_so_quy_uoc_v1_2_20261111.json`.
- Bộ dựng: `bo_dung_chung/ma3/botri_chung_v122u5e.py` (ĐV, U5b), `..._v122u5x.py` (CX), `botri_chung_v123u5c.py` + `botri_cucbo_v123u5b.py` (TV), `..._v122u5n14.py` (N14).
- Đặc tả: `dac_ta/v12/dong_vat_v12.json`, `dong_vat_v12_khe0.8.json`, `thuc_vat_v12.json`.

## 0. Luật bộ kiểm v3 chấm (sự thật v1.2)

Động vật (20): 1.1, 1.2a, 1.3, 1.4a, 1.5, 4.2a, 4.3a, 4.3b, 4.4a, 4.5a, 4.5b, 4.8b, 4.6a, 4.6b, 4.8a, M.a, M.d, M.e, M.c₀, M.h.
Thực vật (28): 2.1, 2.2a, 2.3a, 2.4a, 2.4b, 2.6b, 2.7b, 2.7c, 2.7d, M.d, M.h, 5.1a, 5.2a, 5.3e, 5.4a, 5.6a, 5.7a, 5.9, 10.1a, 10.2a, 10.3a, M.a, 2.7f, 2.6e, 10.2b, 10.2c, 10.3b, M.c₀.
(kiem3.py:620 assert tập id = tập luật trong tệp v1.2.) Ngoài ra `kiem_tham_so_cu` chấm các luật quy ước cũ v4.1/v1.4 chỉ ở phán quyết tham số.

Hằng chung: ε = 1e-6 µm³, chung.tiep_xuc_toi_da_um = 0,01 µm, dung_sai_khoang = 0.

Mục 2 = bảng ĐV (U5b, CX = U5b), mục 3 = khác biệt CX/N14, mục 4 = bảng TV, mục 5 = điều chưa công bố.

## 1. Ghi chép bộ dựng động vật U5b (`botri_chung_v122u5e.py`, đặc tả `dong_vat_v12.json`) — bản nháp

Cơ chế (mã): luật CỨNG = hàng chiếu PROJECT (chứa: dòng 1265–1268; rời khác loại/cùng loại: 1273–1293) + phạt cứng trong J (`phat_lon` 405–425, `phat_nho` 434–449); luật MỀM = `phat_mem_lon` 450–461 (lech_tam, gan), `phat_mem_nho` 482–499 (gan_nho, huong); tầng 2 giải tích `dung_tang2` 1187 (cờ `--tang2`). Lưới h1 = 2·11,5/128 = 0,1797 µm (không có khoá `luoi` trong đặc tả ⇒ toàn phần, dòng 109–125). Mesh = marching cubes từng trường ψ (vobao 1416–1428).
Đặc tả: chứa nhân–màng T=1,0; bào quan–màng T=0,5; nhân con–nhân T=0,3; rời khác loại S_roi=0,3 (khe0.8: 0,8 — khác biệt DUY NHẤT của `dong_vat_v12_khe0.8.json`); cùng loại: ti thể 0,3, lysosome 0, nhân con 0; ER = vành 0 ≤ φ_nhân ≤ D_er=2 µm trừ nón khuyết hướng trung thể (cos 0,55) (`phi_vanh` 367–369); ER–nhân thuộc `mien_cap` ⇒ chỉ hàng u_a+u_b ≥ τ=1e-6, KHÔNG cộng β (1278–1281).
Không có trong đặc tả ĐV: bộ phận `mang_nhan`, `lo_nhan`, `mao`, `ribosome` (grep đặc tả: 0 lần) ⇒ M.d, M.e không dựng.

Bằng chứng vận hành (quan sát, CHƯA rà): (a) 15 ô `bo_kiem_v3/cham_ra/vstdgu5_mu_sau__*/kiem3.json.gz` (U5, mù): 4.5b TRƯỢT 15/15, 4.6b TRƯỢT 15/15, 4.3b TRƯỢT 14/15, 4.5a & 4.8a & M.d & M.e KHÔNG KIỂM ĐƯỢC 15/15, 4.2a & 4.4a MỘT PHẦN 15/15; còn lại ĐẠT. Ca mẫu: Golgi khe nhân 3,14 µm > khe màng 1,44 µm; tâm trung thể ngoài bao lồi Golgi; khe Golgi–trung thể 1,112 > 1,0 (trượt cả tham số); trung vị u lysosome 0,789; tiết diện ti thể 0,163–0,191 µm. (b) 463 tệp `bo_dung_chung/tang_dg/u5v2/doc_lap/*.json` (mọi nhánh U5b/P4/P8/a1/l0/SA un1x/thử–loại nn1x/td, CX, vùng D1/D4/D4A2/D5/k08): 4.5b trượt 463/463, 4.6b 462/463, 4.3b ở hầu hết nhánh vùng D4/D5/k08, M.c₀ chỉ cặp ER–nhân (49 ô, gần hết ở vùng D4A2), 4.8b 4 ô, 4.2a/M.h 1 ô (a1l3).

## 2. Bảng ĐỘNG VẬT — U5b (mọi nhánh P4/P8/thử–loại/SA dùng CÙNG đặc tả + hàng chiếu; khác bộ tối ưu)

Ký hiệu: EXACT = cùng vị từ (hoặc mạnh hơn) được dựng/ràng buộc; EXACT-TS = đúng theo dựng (dải tham số ⊂ dải sự thật), không có ràng buộc; SURROGATE = số hạng mềm khác vị từ; NOT BUILT; POST-HOC = chỉ đạt/hỏng qua PROJECT/marching cubes. "kiem3:N" = dòng kiem3.py; "u5e:N" = dòng botri_chung_v122u5e.py.

| Luật | Vị từ bộ kiểm (sự thật) | U5b | Ghi chú / ngưỡng |
|---|---|---|---|
| 1.1 | 1 màng, 1 nhân, chất TB > ε, nhân trong màng (kiem3:193–200) | EXACT-TS | vỏ `cau_nhieu` 1 cái, nhân `so_luong`=1 (đặc tả `bo_phan[0]`) |
| 1.2a | mọi cạnh OBB màng ≥ 10 µm (TS: ≤ 30) | EXACT-TS | R_vo=10, nhiễu 0,25 ⇒ ≈20 µm; không ràng buộc |
| 1.3 | 1 nhân khối đặc liền trong TB | EXACT | chứa nhân–màng T=1,0 (quan_he.chua[0]; hàng u5e:1265–1268); lệch tâm = mềm (u5e:453–454) |
| 1.4a | nhân − màng ≤ ε | EXACT (mạnh hơn) | T_nhan_mang = 1,0 µm + β |
| 1.5 | tỉ số Ø nhân/TB 0,35–0,6; mọi cạnh OBB nhân 5–10 µm | EXACT-TS | dai_nhan 3,9–5,0 × trục [1; 0,943; 0,896] ⇒ cạnh 6,99–10,0; tỉ số ≈ 0,37–0,47. Cận trên 10,0 = đúng biên sự thật, DS = 0 ⇒ rủi ro POST-HOC do lưới |
| 4.2a | 1–5 nhân con, mỗi cái 0,5–5 µm, trong khối nhân | EXACT (+EXACT-TS cỡ) | so_nhan_con integers(1,4) = 1–3; Ø 0,8–1,8; chứa trong nhân T_nhancon=0,3. "Không màng bao" → MỘT PHẦN |
| 4.3a | ti thể ⊂ màng, không giao nhân | EXACT | chứa T_mang 0,5; roi_vung_cam (ER∪nhân, s=0,3); roi_chieu nhân–ti thể |
| 4.3b | độ dày cục bộ 0,2–1 µm | EXACT-TS trên khuôn, POST-HOC trên mesh | khuôn ống rc ∈ {0,25…0,45}×nửa dài 0,5–1,0 ⇒ Ø 0,25–0,9 trên giấy; Ø < 2h1 = 0,36 µm là dưới lưới ⇒ marching cubes làm mỏng còn 0,16–0,19 → TRƯỢT. Không có hàng/phạt nào cho độ dày |
| 4.4a | ER chạm mặt ngoài nhân: khe ≤ 0,01, diện tích chạm > 0 | EXACT (qua cấu tạo) | vành `er_trong`=0 ⇒ mặt trong ER ≡ mặt nhân (u5e:367–369); hàng ER–nhân τ=1e-6 không β (u5e:1278–1281). ER hạt/trơn, ribosome, "trải khắp" → MỘT PHẦN (ER chỉ là vành dày 2 µm quanh nhân) |
| 4.5a | Golgi ≥ 2 khối: 4–11 túi, mỗi túi 0,7–1,1 µm | EXACT ở tầng 2, nhưng bộ kiểm KHÔNG thấy | `tang2` chong_dia n∈4–8, Ø 0,7–1,1 (u5e:1205–1215); túi là bộ phận loại `tui_golgi` riêng, Golgi cha vẫn 1 khối elip ⇒ kiem3:280–287 KHÔNG KIỂM ĐƯỢC 15/15 |
| 4.5b | (i) khe Golgi–nhân < khe Golgi–màng; (ii) tâm trung thể trong bao lồi Golgi | SURROGATE (i: KHÔNG có số hạng nào), (ii): KHÔNG KHẢ THI theo cấu tạo | chỉ có `gan` trung thể→Golgi: relu²(φ_golgi(tâm tt) − r_tt − 1,0) (u5e:455–460, dv.gan_tt_golgi=1,0, mềm). Golgi elip lồi + roi_lon S_roi ⇒ trung thể RỜI Golgi ⇒ tâm không bao giờ trong bao lồi. Không có số hạng Golgi–nhân |
| 4.6a | mọi cạnh OBB lysosome 0,1–1,2 | EXACT-TS | cầu r 0,2–0,4 ⇒ Ø 0,4–0,8 |
| 4.6b | trung vị u < 0,5 và ∃ u ≥ 0,5 | SURROGATE (yếu, ngược chiều trên thực tế) | gan_nho: 0,05·relu²(φ_nhân(TÂM lysosome) − 4,0) (u5e:484–486; khoảng TÂM–mặt, không mặt–mặt); không số hạng cho ∃u ≥ 0,5. Vành ER 2 µm + khe 0,3 cấm vùng gần nhân ⇒ trung vị u ≈ 0,79 |
| 4.8a | 2 trung tử 0,3–0,5 × 0,15–0,25, góc ≥ 75° | EXACT ở tầng 2, bộ kiểm KHÔNG thấy | `hai_tru_vuong_goc` 90° theo cấu tạo (u5e:1194–1202); trung thể cha 1 khối cầu ⇒ kiem3:336–337 KKD 15/15 |
| 4.8b | khe trung thể–nhân < khe trung thể–màng | SURROGATE | gan: relu²(φ_nhân(tâm) − r − 1,5) (dv.gan_tt_nhan=1,5) + nón khuyết ER hướng trung thể; đạt 15/15 nhưng trượt 4/463 |
| M.a | mọi bộ phận − màng ≤ ε | EXACT (mạnh hơn) | T_mang 0,5 |
| M.d | vỏ nhân hai lớp + ≥ 1 lỗ | NOT BUILT | không có `mang_nhan`/`lo_nhan` trong đặc tả ĐV (dv.mang_nhan_ve có trong tệp tham số nhưng không dùng) ⇒ KKD |
| M.e | ≥ 1 mào mỗi ti thể | NOT BUILT | không có `mao` ⇒ KKD |
| M.c₀ | mọi cặp khác loại: độ sâu xuyên ≤ 0,01 | EXACT (mạnh hơn, s=0,3/0,8 + β) trừ ER–nhân: SURROGATE lưới (τ=1e-6, không β) ⇒ POST-HOC | ER–nhân là cặp trượt duy nhất (49/463, chủ yếu D4A2) |
| M.h | cùng loại: giao ≤ ε | EXACT (mạnh hơn) | ti thể 0,3; lysosome 0; nhân con 0 (+β) (u5e:1283–1293) |

Thay đổi theo nhánh: `--khong_mem` bỏ mọi số hạng mềm (4.5b/4.6b/4.8b/1.3 lệch tâm) → không còn cả surrogate. SA (`--cach u_nhiet`) và thử–loại (`--cach ngau_nhien`) dùng cùng J/S(g) (u5e:53–64), cùng hàng PROJECT. `khe0.8` chỉ đổi S_roi 0,3 → 0,8 (đặc tả dòng 1203).

## 3. CX (`botri_chung_v122u5x.py`) và N14 (`botri_chung_v122u5n14.py`) — khác U5b ở đâu

**CX** = u5e + `--cx` (PROJECT xấp xỉ trong vòng Adam + neo tư thế, u5x:1221–1315; u5x:105 cấm dùng chung `--u5`). Không đổi đặc tả, hàng chiếu hay số hạng mềm ⇒ cột CX = cột U5b ở MỌI luật. Kiểm chứng (quan sát, chưa rà): 30 ô `vstdgu5x*` — 4.5b, 4.6b TRƯỢT 30/30; 4.3b TRƯỢT 29/30; 4.5a/4.8a/M.d/M.e KKD 30/30 (giống U5b `vstdgu5r1`: 4.3b 19/30, M.c₀ 5/30).

**N14** = u5f + vật chứa THẬT (Allen hiPSC): đặc tả `bo_dung_chung/n14/dac_ta/dong_vat_n14_<CellId>.json` (5 tế bào; `n14d5/dac_ta/dong_vat_n14k08_*` cho vùng k08), sinh từ `dac_ta/v12/dong_vat_v12.json` (md5 8d7f6bfe). Khác biệt (khoá `n14_thay_doi`): vỏ `luoi` (SDF tới mesh màng thật), nhân `co_dinh` (mesh nhân thật, tư thế + cỡ cố định, dai_nhan = 6,77 µm), `lech_tam` = [] (bỏ), khởi tạo/đặt lại `vung_mau`. `vanh`, `tang2`, `chiem` giữ nguyên; luật + tham số khe giữ nguyên. Hệ quả cho bảng:
- 1.2a, 1.5, 1.3 (lệch tâm), 1.4a của nhân: KHÔNG do bộ dựng quyết — là DỮ LIỆU VÀO (Allen). 90 ô `vsn14chinh`: 1.2a TRƯỢT 18/90, 1.5 TRƯỢT 45/90 — bộ dựng không có cơ chế sửa.
- M.c₀ có nhân: nhân thật có thể không genus 0 ⇒ theo kiem3:107–147 (§9) cặp có nhân thành KKD; 19/90 ô M.c₀ TRƯỢT.
- 4.8b TRƯỢT 10/90 (khác ĐV tổng hợp 4/463): surrogate gan_tt_nhan 1,5 tính trên tâm trung thể tới nhân thật; GIẢ THUYẾT chưa kiểm: màng thật dẹt (nửa bao 14,4 × 7,5 × 5,7 µm) ⇒ khe màng nhỏ hơn.
- 4.3b TRƯỢT 10/90 (dưới lưới như U5b); 4.5b 90/90, 4.6b 89/90 trượt; 4.5a/4.8a/M.d/M.e KKD 90/90 — như U5b.

## 4. Bảng THỰC VẬT — `botri_chung_v123u5c.py` + `botri_cucbo_v123u5b.py`, đặc tả `thuc_vat_v12.json`

Đặc tả TV: `pham_vi` = "tầng 1: vách, màng, nhân, lục lạp, ti thể, Golgi, peroxisome, lưới nội chất, không bào"; KHÔNG có khoá `tang2` (grep v123u5c: chỉ mã tầng 2 chung, không mục nào cho TV). Lưới cục bộ h = 0,05 µm (`luoi`). Vỏ `hop_bo_tron` nửa cạnh 40×16×16, bo 3, lớp vách 1,0–1,5 + màng 0,15–0,25. Chứa trong màng: nhân 0,3; lục lạp 0,3; ti thể/Golgi/peroxisome 0,3. Rời khác loại 0,3; cùng loại: lục lạp/ti thể/peroxisome 0,3, Golgi 2,0. Mềm: `gan` nhân–màng ≤ 0,5 (tâm−mặt trừ bán kính, trọng số 1); `gan_nho` (đo bằng mẫu mặt, w_gan_cung = 100) lục lạp–màng ≤ 0,5, ti thể/Golgi/peroxisome–(màng hoặc nhân) ≤ 1,0; `huong` lục lạp lệch ≤ 30°, w_huong = 10 (v123u5c:502–503, trục z khuôn · pháp tuyến vỏ TẠI TÂM). ER: vùng phụ thuộc quanh nhân `trong`=0, dày 0,5–1,0, ống ngoại vi r 0,1, trừ mọi bộ phận + 0,3, chạm tối thiểu 1,0 µm². Không bào: vùng phụ thuộc, tỉ lệ thể tích 0,68–0,85, liền khối, cách màng 1,0, cách bộ phận 0,3. `ngan_sach.giam = [luc_lap]`: số lục lạp bị GIẢM khi ước lượng chỗ trống < ngưỡng (cùng mã u5e:985–990).

Bằng chứng (quan sát, CHƯA rà): 25 ô `bo_kiem_v3/cham_ra/vsu5tv*_mu_sau__*` (lớp Sinh10). Trượt: 5.1a, 5.2a, 5.3e, 10.1a, 10.2a 25/25; 2.7f 15/25 (418 ca, 49–54°); 2.7c 16/25 (281 ca, độ lõm 0,0221 > 0,02); M.d 17/25 (255 ca, khe 1,10 > 1,0); 2.3a 15/25 (giao màng–vách 1,4e-5 > ε 1e-6 µm³); 2.7b 13/25 (71 ca, khe 0,575–0,589 > 0,5). KKD: 10.2b, 10.2c, 10.3b 25/25. 10.3a ĐẠT với 0 granum (đạt rỗng).

| Luật | Vị từ bộ kiểm (sự thật) | Bộ dựng TV | Ghi chú / ngưỡng |
|---|---|---|---|
| 2.1 | 1 thành, 1 màng, ≥1 nhân/không bào/lục lạp | EXACT-TS | vỏ 2 lớp, nhân 1, KB 1, lục lạp integers(30,51) (có thể bị `ngan_sach` giảm, vẫn ≥ 30) |
| 2.2a | V trong màng 60–100 pL | EXACT-TS (gần cận dưới) | hộp 80×32×32 − 2(vách+màng) ⇒ ≈ 62–68,5 pL (đo 67,08); không ràng buộc |
| 2.3a | màng trong thành, giao vật liệu thành ≤ ε | EXACT trên trường, POST-HOC trên mesh | khe màng–vách = 0 (tv.khe_mang_vach = 0) ⇒ hai mesh chung mặt ⇒ giao 1,4e-5 µm³ > ε ⇒ trượt 15/25 khi không sửa lưới |
| 2.4a | đúng 1 KB, 1 khối | EXACT | `lien_khoi`: true; "màng riêng" → MỘT PHẦN |
| 2.4b | V(KB)/V(trong màng) 0,68–0,85 | EXACT | `ti_le_the_tich` [0,68; 0,85] (vung_phu_thuoc[1]) |
| 2.6b | khe nhân–màng ≤ 0,5; không trong KB | SURROGATE (cùng ngưỡng 0,5, mềm) + EXACT cận dưới 0,3 | `gan` tâm−bán kính (đúng cho cầu); KB trừ mọi bộ phận ⇒ không trong KB: EXACT |
| 2.6e | V nhân / V trong màng 0,0016–0,023 | EXACT-TS (sát cận dưới) | nhân cầu r 3,0–3,5 ⇒ 0,0016–0,0029; r = 3,0 với vách mỏng ≈ 0,00165 |
| 2.7b | mỗi lục lạp khe ≤ 0,5; rải ≥ 3 mặt hộp | SURROGATE (khe: dải 0,3 cứng – 0,5 mềm w=100); "≥ 3 mặt": NOT BUILT (không số hạng) | trượt 13/25 |
| 2.7c | lồi ≥ 0,98, lõm ≤ 0,02 µm, tròn ≥ 0,9, hồ sơ 0,48–0,6 | EXACT trên khuôn `thau_kinh`, POST-HOC trên mesh | lõm do lưới 0,05 vượt 0,02 ⇒ trượt 16/25; KHÔNG được đánh dấu surrogate trong bài |
| 2.7d | Ø 4–10, dày 1,11–2,5 | EXACT-TS | khuôn thau_kinh0–2: Ø 5–7, dày 2,0–2,5 (thau_kinh1/2 chạm đúng 2,5 = cận sự thật, DS = 0) |
| 2.7f | trục ngắn PCA lệch pháp tuyến mặt màng gần nhất ≤ 30° | SURROGATE | `huong` w=10: trục z khuôn vs pháp tuyến vỏ tại TÂM (sai phân φ_vo); bộ kiểm: PCA + pháp tuyến mặt gần nhất; trượt 15/25 |
| M.d | ti thể/Golgi/peroxisome cách màng hoặc nhân ≤ 1,0 | SURROGATE (cùng ngưỡng, mềm w=100, đo mẫu mặt) | trượt 17/25 |
| M.h | cùng loại giao ≤ ε | EXACT (mạnh hơn) | 0,3 / Golgi 2,0 |
| M.a | mọi bộ phận trong khối trong màng | EXACT | chứa 0,3 |
| M.c₀ | khác loại: độ sâu xuyên ≤ 0,01 | EXACT (s 0,3 + β); ER/KB: EXACT qua phép trừ `tru: *` + cách 0,3 (ER–nhân: `mien_cach` nhân, `trong`=0 ⇒ chạm) | đạt 25/25 |
| 5.1a | ≥ 1 nhân con trong nhân | NOT BUILT | không có `nhan_con` ⇒ TRƯỢT 25/25 (không phải KKD) |
| 5.2a | vỏ nhân 2 lớp, ≥ 1 lỗ | NOT BUILT | TRƯỢT 25/25 |
| 5.3e | ≥ 1 mào mỗi ti thể | NOT BUILT | TRƯỢT 25/25 (2 037 ca) |
| 5.4a | ER chạm mặt ngoài màng nhân, diện tích > 0 | EXACT | ER `trong`=0 quanh nhân (không có màng nhân ⇒ bộ kiểm dùng mặt nhân); ER hạt/trơn → MỘT PHẦN |
| 5.6a | ≥ 2 Golgi, không cặp nào chạm | EXACT | 5–10 thể, khe cùng loại 2,0 |
| 5.7a | mọi cạnh OBB peroxisome 0,1–1 | EXACT-TS | cầu r 0,25–0,5 ⇒ Ø 0,5–1,0 (cận trên = cận sự thật, DS = 0) |
| 5.9 | 0 trung thể | EXACT | không có loại trung_the |
| 10.1a | 1 màng trong mỗi lục lạp | NOT BUILT | TRƯỢT 25/25 |
| 10.2a | ≥ 1 granum mỗi lục lạp | NOT BUILT | TRƯỢT 25/25 |
| 10.2b / 10.2c / 10.3b | cỡ granum / số thylakoid / tấm nối | NOT BUILT | KKD 25/25 |
| 10.3a | đồ thị granum–tấm nối liền | NOT BUILT | ĐẠT RỖNG (0 granum ⇒ "liền") |

`botri_cucbo_v123u5b.py` chỉ là phần chiếu/mesh lưới cục bộ (+ `--giu_truoc`, đo HD); không thêm vị từ luật.

## 5. Điều bất ngờ CHƯA nêu trong `main.tex` (Method §surrogate dòng 152–156, §checker 240–248, Limitations 431–440)

Bài đã nêu: surrogate cho 4.5b/4.6b/4.8b/2.7b/2.7f/M.d với ngưỡng tuyệt đối; vế "cạnh nhân" của 4.5b không có số hạng; sửa lưới cho cặp chạm (ER–nhân, màng–vách); "một số luật mô tả cấu trúc dưới bước lưới và không được dựng". Chưa nêu:

1. **4.5b vế (ii) KHÔNG KHẢ THI theo cấu tạo, không chỉ là surrogate kém.** Golgi là elip lồi (khuôn `golgi`, `ai_bo_golgi` bỏ khuôn AI) và trung thể bị ép RỜI Golgi (roi_lon S_roi 0,3/0,8 + β) ⇒ tâm trung thể không bao giờ nằm trong bao lồi Golgi ⇒ 4.5b trượt 463/463 ô ĐV, 90/90 N14, 30/30 CX bất kể phương pháp. Bài viết "centrosome–Golgi distance ≤ 1,0 for the second clause" — một số hạng khoảng cách không thể thoả "bao quanh". Ngoài ra số hạng dùng khoảng TÂM trung thể→mặt Golgi trừ cỡ (u5e:459–460), còn phán quyết tham số của bộ kiểm dùng khe mặt–mặt chính xác; ca đo 1,112 > 1,0 trượt cả tham số.
2. **4.6b ngược chiều trên thực tế:** trung vị u ≈ 0,79 (lysosome nằm phía màng). Số hạng = 0,05·relu²(φ_nhân(TÂM) − 4,0), khoảng TÂM–mặt (không trừ bán kính) trong khi tham số bộ kiểm là trung vị khe mặt–mặt ≤ 4,0; vành ER 2 µm + khe 0,3 cấm vùng sát nhân; vế "∃ u ≥ 0,5" không có số hạng. Trượt 462/463.
3. **4.5a và 4.8a dựng ĐÚNG ở tầng 2 nhưng bộ kiểm KHÔNG chấm được:** túi Golgi (4–8, Ø 0,7–1,1) và hai trung tử vuông góc là bộ phận con riêng (`tui_golgi`, `trung_tu`); bộ kiểm đếm khối trên mesh CHA (một khối) ⇒ KHÔNG KIỂM ĐƯỢC ở 100 % ô. Bài không nói hai luật này luôn KKD.
4. **Luật không dựng không chỉ vì "dưới bước lưới":** ĐV thiếu vỏ nhân/lỗ nhân (M.d) và mào (M.e) ⇒ KKD 100 %. TV chỉ dựng tầng 1 ⇒ nhân con (Ø 1–2 µm, không dưới lưới h = 0,05), vỏ nhân, mào, màng trong lục lạp, granum KHÔNG dựng, và vì chấm lớp Sinh10 nên 5.1a, 5.2a, 5.3e, 10.1a, 10.2a là TRƯỢT 25/25 (không phải KKD); 10.2b/10.2c/10.3b KKD; 10.3a ĐẠT RỖNG (0 granum). Cần ghi rõ phạm vi tầng 1 hoặc chấm TV ở KHTN6 / loại các luật này khỏi bảng.
5. **4.3b (tiết diện ti thể) trượt do lưới, không do phương pháp:** khuôn ống cho Ø 0,25–0,9 µm trên giấy, nhưng Ø < 2h1 = 0,36 µm là dưới lưới; mesh đo 0,16–0,19 µm < 0,2 ⇒ trượt 14/15 (U5), hầu hết ô D4/D5/k08; 0 ô ở D1/D4A2. Bộ phân tích còn có trường `V_test_tru_43b_duoi_luoi` (loại 4.3b dưới lưới khỏi V_test) — nếu dùng thì phải công bố.
6. **M.c₀ ER–nhân:** cặp duy nhất trượt M.c₀ (49/463 ô, chủ yếu D4A2); hàng chiếu của cặp này là u_a+u_b ≥ τ = 1e-6 KHÔNG cộng β (u5e:1278–1281) trong khi mọi cặp khác có s + β ⇒ bất cân xứng có chủ đích nhưng chưa nêu; bảng W1 lại chấm KHÔNG sửa lưới.
7. **Trượt do lưới ở TV chưa được đánh dấu:** 2.3a (khe màng–vách = 0 ở bộ dựng vs ε = 1e-6 µm³ ở bộ kiểm; 15/25), 2.7c (độ lõm 0,0221 > TV_NEP_SAU 0,02; 16/25) — 2.7c không mang dấu † nên dễ bị đọc là lỗi phương pháp.
8. **2.7b vế "rải trên ≥ 3 mặt" không có số hạng**; 2.7f surrogate đo pháp tuyến vỏ TẠI TÂM lục lạp + trục z khuôn, bộ kiểm đo PCA + pháp tuyến mặt màng gần nhất (định nghĩa khác, không chỉ ngưỡng khác); trọng số 10 so với 100 của luật gần.
9. **Dải tham số chạm đúng cận sự thật với DS = 0** (rủi ro trượt do lưới): 1.5 cạnh nhân tới 10,0; 5.7a Ø peroxisome tới 1,0; 2.7d dày lục lạp tới 2,5; 2.6e sát cận dưới 0,0016; 2.2a ≈ 62–68,5 pL sát 60.
10. **N14: 1.2a (18/90) và 1.5 (45/90) trượt vì dữ liệu Allen**, không do bộ dựng (vỏ + nhân cố định, `lech_tam` bị bỏ); cần tách khỏi so sánh phương pháp. Nhân thật có thể không genus 0 ⇒ M.c₀ cặp có nhân thành KKD (kiem3 §9).
11. **Số lục lạp TV thích nghi:** `ngan_sach.giam = [luc_lap]` giảm số lục lạp đã rút khi ước lượng chỗ trống thấp (cùng phép cho mọi nhánh) — chưa nêu ở Configurations.
12. `--khong_mem` (đối chứng) bỏ cả surrogate ⇒ nếu có bảng dùng nhánh này, các dòng † không còn số hạng nào.

## 6. Trạng thái
CHƯA rà soát chéo. Mọi tỉ lệ trượt/KKD ở trên là đếm trực tiếp từ `kiem3.json.gz`/`kiem3.csv` (`bo_kiem_v3/cham_ra/`) và `u5v2/doc_lap/*.json`, chưa đối chiếu với bảng trong bài; giải thích nguyên nhân (mục 2, 5, 10, N14 4.8b) là suy từ mã, cần agent độc lập kiểm lại trước khi đưa vào bài.
