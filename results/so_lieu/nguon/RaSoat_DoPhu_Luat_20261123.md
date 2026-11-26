# Rà soát chéo độc lập: DoPhu_Luat_BoDung_20261123.md (23/11/2026)

Trạng thái: XONG một vòng rà độc lập. Agent rà chỉ đọc server và tự đếm lại từ dữ liệu gốc (`bo_kiem_v3/cham_ra/*/kiem3.json.gz`, `loai.json`, `mesh.npz`, `bo_dung_chung/tang_dg/u5v2/doc_lap/*.json`). Mã đọc trực tiếp: `bo_kiem_v3/kiem3.py`, `bo_dung_chung/ma3/botri_chung_v122u5e.py` (u5e), `v122u5x`, `v122u5n14`, `v123u5c`. Đặc tả: `dac_ta/v12/dong_vat_v12.json`, `thuc_vat_v12.json`. Tài liệu được rà KHÔNG được dùng làm tiền đề.

## Tóm tắt phán quyết

| # | Mệnh đề | Phán quyết | Ảnh hưởng so sánh phương pháp? |
|---|---|---|---|
| 1 | 4.5b vế "bao quanh" không bao giờ thoả | XÁC NHẬN | Không. Mọi nhánh đều trượt, chỉ hạ tỉ lệ đạt tuyệt đối |
| 2 | 4.6b trung vị u ≈ 0,79 | XÁC NHẬN | Không đáng kể (trượt 462/463) |
| 3 | 4.5a, 4.8a dựng ở tầng 2 nhưng bộ kiểm luôn KKD | XÁC NHẬN | Không (KKD ở mọi ô) |
| 4 | TV chỉ dựng tầng 1: 5 luật trượt 25/25, 10.3a đạt rỗng | XÁC NHẬN | Không (mọi nhánh TV); chỉ ảnh hưởng tuyệt đối |
| 5 | Do lưới: 4.3b, TV 2.3a, 2.7c | MỘT PHẦN | 4.3b CÓ (tần suất khác nhau theo nhánh, và nó nằm trong V_test); 2.7c chưa rõ |
| 6 | M.c₀ ER–nhân: không cộng β, lỗi bộ dựng? | MỘT PHẦN (không cộng β là đúng; gọi là lỗi là SAI) | Gần như không. Cùng mã ở mọi nhánh, độ sâu xuyên do nhân quyết định |
| 7 | N14: 1.5 trượt 45/90, 1.2a trượt 18/90 do dữ liệu Allen | XÁC NHẬN, kèm sắc thái quan trọng | Không giữa các phương pháp; CÓ nếu so vật chứa thật với vật chứa tổng hợp |

## 1. 4.5b: XÁC NHẬN
- Bộ kiểm (kiem3.py:289–307) có hai vế. Vế (i): khe Golgi–nhân < khe Golgi–màng. Vế (ii): tâm khối trung thể nằm trong bao lồi (ConvexHull + Delaunay) của các đỉnh mesh Golgi CHA.
- Đếm lại 15 ô `vstdgu5_mu_sau`: vế (ii) sai 15/15, vế (i) cũng sai 15/15. Ca mẫu: khe nhân 3,138 so với khe màng 1,443, khe Golgi–trung thể 1,112. Trong 463 tệp doc_lap: 4.5b trượt 463/463. N14: 90/90.
- Vì sao không thoả được: khuôn `golgi` là `elip` (đặc tả `khuon`). Khuôn AI bỏ Golgi (`ai_bo_golgi = ['golgi']`), nên Golgi luôn lồi và bao lồi gần trùng chính khối Golgi. Trung thể (cầu, `khuon_trung_the = cau`) cùng Golgi nằm trong `roi_lon` và `roi_chieu`, nên hàng PROJECT buộc chúng cách nhau ≥ s + β (u5e:1280–1285; S_roi 0,3, β = √3·h1 = 0,311). Muốn tâm trung thể nằm trong bao lồi thì tâm đó phải nằm trong khối Golgi, tức xuyên ≥ r_tt ≈ 0,42 µm. Hàng chiếu cấm điều đó. Vậy vế (ii) không khả thi theo cấu tạo ở MỌI nhánh. Túi tầng 2 nằm trong Golgi nên cũng không đổi được bao lồi.
- Số hạng mềm duy nhất (u5e:460, `gan` trung thể→Golgi ≤ 1,0, tính từ tâm trừ cỡ) không thể tạo quan hệ "bao quanh". Vế (i) không có số hạng nào, nhưng nó KHÔNG bất khả thi theo cấu tạo: Golgi có thể vào nón khuyết của ER.

## 2. 4.6b: XÁC NHẬN
- 15 ô U5: trung vị u từng ô nằm trong 0,747–0,796, trung vị chung 0,789. Vế "∃ u ≥ 0,5" đạt (ví dụ 0,876).
- Số hạng: u5e:486, `0,05·relu²(φ_nhân(TÂM) − 4,0)`. Khoảng cách tính từ tâm, không trừ bán kính, vì luật không có khoá `do`. Không có số hạng cho vế ∃. Lysosome thuộc `roi_vung_cam` của ER (s = 0,3), nên mặt lysosome cách nhân ≥ D_er + 0,3 = 2,3 µm (ngoài nón khuyết). Phần lớn vùng cho phép vì vậy có u ≈ 0,45 hoặc lớn hơn. Tham số bộ kiểm "trung vị khe ≤ 4" cũng trượt (4,283).
- 4.6b trượt 462/463. Ô đạt duy nhất thuộc e3b_B_D1_un1x. Mọi nhánh đều trượt, nên không ảnh hưởng so sánh.

## 3. 4.5a, 4.8a: XÁC NHẬN
- Bộ dựng: `dung_tang2` (u5e:1187–1232) tạo các mesh con riêng mang loại `tui_golgi` (chong_dia, 4–8 túi) và `trung_tu` (hai_tru_vuong_goc). Chúng được ghi vào `mesh_vo_bao.npz` (u5e:1433–1438). Ô mẫu `vstdgu5_mu_sau__0b6ee6a83fed`: mesh và `loai.json` có tui_golgi 8, trung_tu 2, nghĩa là bộ kiểm CÓ nhận các bộ phận này.
- Bộ kiểm: 4.5a (kiem3.py:275–288) chỉ duyệt `theo["golgi"]` và đếm `b.so_khoi` của mesh CHA. Một khối thì gán `kkd`. 4.8a (kiem3.py:336–353) chỉ duyệt `theo["trung_the"]`; `so_khoi == 1` thì KKD. Trong kiem3.py, `tui_golgi`/`trung_tu` chỉ xuất hiện ở `TO_TIEN_DV` (dòng 152, dùng để loại cặp tổ tiên–con ở M.c₀). Không có `con_loai` cho động vật. Thực vật thì có (`tb.con_loai`, dòng 551, cho 10.1a). Đây là khe hở của bộ kiểm, không phải của bộ dựng.
- 15/15 ô U5 đều KKD ở cả hai luật. Phụ: tham số 4.8a trượt (bán kính đo 0,411 < 0,42), nhiều khả năng do marching cubes làm co.

## 4. Thực vật chỉ tầng 1: XÁC NHẬN
- `thuc_vat_v12.json`: `pham_vi` = "tầng 1: …", không có khoá `tang2`. Bộ phận chỉ gồm nhân, lục lạp, ti thể, Golgi, peroxisome, ER và không bào. Không có nhan_con, mang_nhan/lo_nhan, mao, mang_trong_lap, granum hay tấm nối. Thêm nữa, với lưới cục bộ, u5e:1425–1426 `exec` botri_cucbo rồi `sys.exit`, nên `--tang2` cũng không chạy.
- Đếm lại 25 ô `vsu5tv*_mu_sau` (lớp Sinh10): 5.1a, 5.2a, 5.3e, 10.1a, 10.2a TRƯỢT 25/25. 10.2b, 10.2c, 10.3b KKD 25/25. 10.3a ĐẠT 25/25 với `gia_tri = 0` thành phần (đạt rỗng).
- Phần còn thiếu: 5.1a thiếu nhân con trong nhân (đếm 0 < 1). 5.2a: "mang_nhan khong co → phai co" (vỏ nhân hai lớp và lỗ nhân). 5.3e: 0 mào ở mỗi ti thể. 10.1a: 0 màng trong ở mỗi lục lạp. 10.2a: 0 granum. Tất cả đều là cấu trúc con BÊN TRONG một bào quan (tầng 2). Riêng nhân con (Ø 1–2 µm) không hề dưới lưới h = 0,05, nên lời giải thích "dưới bước lưới" không đúng cho 5.1a.
- Bất cân xứng bộ kiểm: ở động vật, các bộ phận vắng (M.d, M.e) bị chấm KKD; ở thực vật, các bộ phận vắng bị chấm TRƯỢT.

## 5. Các luật trượt "do lưới": MỘT PHẦN
- **4.3b:** cơ chế XÁC NHẬN. `khuon_ong` có rc 0,25–0,45 × nửa dài 0,5–1,0, cho Ø trên giấy 0,25–0,9 µm. h1 = 0,1797, ngưỡng dưới lưới 2·h1 = 0,36. Nhưng dải đo sai: trong 15 ô U5, tiết diện của các ti thể trượt là **0,107–0,191 µm**, không phải 0,16–0,19 (con số đó chỉ là một ô mẫu).
  - Mệnh đề "không do phương pháp" là SAI một phần. Trong cùng vùng, tần suất khác nhau theo nhánh. D5: e3b nn1x 10/10 ô trượt và 7–20 ti thể dưới lưới mỗi ô; un1x 8/10 (0–9 mỗi ô); P4 3/10, l0 3/10, P8 4/10, U5b-r 4/10 (0–2 mỗi ô).
  - Trường `V_test` trong doc_lap CÓ đếm các ca này. `V_test_tru_43b_duoi_luoi` = V_test trừ ca dưới lưới; ở D5, nn1x có 179 so với 50 sau khi trừ. Nếu bảng so sánh dùng V_test thì 4.3b làm lệch có lợi cho PLACE. Phải nói rõ bảng dùng trường nào.
  - D1 và D4A2: 0 ô trượt 4.3b. k08/K4/K8: gần 100 %.
- **TV 2.3a:** XÁC NHẬN là lỗi mesh. Cả 15 ô trượt đều thuộc lô `u5tv/khoi` (không sửa mesh) và cùng một giá trị giao 1,4e-5 µm³. 10 ô lô `cham/lo_1…8` có `sua_mesh.json` (sửa cặp vách|mang) và đều đạt. Hệ quả: khác biệt giữa các nhánh TV ở 2.3a phản ánh chính sách sửa mesh, không phản ánh phương pháp. Cần kiểm mọi nhánh có cùng chính sách.
- **TV 2.7c:** SAI một phần. 281/2 274 lục lạp trượt. Độ lõm trung vị 0,042 µm; 116 ca > h = 0,05; 38 ca > 0,1; tối đa 0,328 µm. Các vết lõm lớn hơn vài bước lưới không thể quy hết cho lưới, có khả năng là móp do phép chiếu hay va chạm. Chưa tách theo nhánh.

## 6. M.c₀ ER–nhân: MỘT PHẦN
- XÁC NHẬN: hàng chiếu của cặp miễn khe (`mien_cap = [[luoi_noi_chat, nhan]]`) có vế phải là `TAU_C0 = tol_con = 1e-6`, KHÔNG cộng β. Vị trí: u5e:1277–1279 (`them_hang(ix[a],1,ix[b],1,0,TAU_C0)`); u5x:1593; u5n14:1393 (`hang_cap(a,b,TAU_C0)`); v123u5c:1413. Mọi cặp khác dùng s + bu (u5e:1283–1285). Mã giống nhau ở mọi bộ dựng và mọi nhánh.
- Gọi đây là "lỗi bộ dựng" là SAI. β = √3·0,1797 = 0,311 µm. Nếu cộng β, ER phải cách nhân ≥ 0,31 µm, phạm 4.4a (khe ≤ 0,01, diện tích chạm > 0). Bỏ β là bắt buộc. Xuyên còn lại đến từ marching cubes của hai trường chạm nhau ở h1 = 0,18.
- Bằng chứng xuyên do lưới và do hình nhân: độ sâu xuyên của 49 ca nằm trong 0,010–0,050 µm (≤ 0,28·h1). 100 % số ca ở vùng D4A2 (khuôn nhân AI) và N14 (nhân thật); vùng có nhân elip (D1/D4/D5/k08) không có ca nào. Cùng seed thì độ sâu trùng nhau giữa các nhánh: 6518 cho 0,02302 ở U5b, P4 và l0. N14: một tế bào cho 0,026756 ở cả 9/9 nhánh, một tế bào khác 0,018457 ở 7 nhánh. Đây là lớp L2 phụ thuộc nhân, không phụ thuộc phương pháp. Số ô trượt theo nhánh ở D4A2 dao động 2–9/10.
- Tài liệu viết "chủ yếu D4A2"; thực tế là 49/49 ở D4A2.
- Tài liệu (mục 3) cho rằng ở N14 "nhân thật không genus 0 nên cặp có nhân thành KKD". Điều này không xảy ra trong 90 ô `vsn14chinh`: không ô nào có ghi chú KKD ở M.c₀. Cả 19 ô trượt đều là cặp ER–nhân, đều trên nhân thật. ⇒ RÚT LẠI mệnh đề này.

## 7. N14: XÁC NHẬN, kèm sắc thái
- 90 ô `vsn14chinh` = 45 ô vật chứa THẬT (5 tế bào Allen 593409, 707681, 798836, 989830, 997092, mỗi tế bào 9 ô) + 45 ô vật chứa TỔNG HỢP (khuôn `nhan`, màng ≈ 20 µm). Tài liệu không nêu điều này.
- 1.5 trượt đúng 45/45 ô thật và 0/45 ô tổng hợp. Tỉ số Ø nhân/TB 0,639–0,677 > 0,6, cạnh OBB nhân tới 19,2 > 10. Giá trị giống hệt nhau ở mọi nhánh của cùng tế bào, nên hoàn toàn do dữ liệu.
- 1.2a trượt 18 = 2 tế bào (593409, 997092) × 9. Cạnh OBB nhỏ nhất của màng là 7,04 và 8,13 < 10 µm. Cũng do dữ liệu.
- M.c₀ 19/90 đều ở ô thật (19/45). 4.8b 10/90: 9 ở ô thật, 1 ở ô tổng hợp. Điều này hợp với giả thuyết màng dẹt, nhưng chưa được kiểm trực tiếp.
- Ảnh hưởng: không ảnh hưởng so sánh giữa các phương pháp trong cùng vật chứa. CÓ ảnh hưởng nếu bài so vật chứa thật với vật chứa tổng hợp, hoặc gộp tỉ lệ đạt.

## Cần sửa trong DoPhu_Luat_BoDung
1. 4.3b: thay dải đo thành 0,107–0,191. Bỏ câu "không do phương pháp". Nêu rõ V_test có hay không trừ 4.3b dưới lưới.
2. M.c₀: ghi "100 % D4A2". Ghi rằng bỏ β là bắt buộc vì 4.4a (không phải lỗi bộ dựng). Rút lại mệnh đề KKD genus ở N14.
3. 2.7c: không gọi là lỗi do lưới thuần.
4. 2.3a: nêu rõ phụ thuộc bước sửa mesh (lô khoi so với lô lo_).
5. N14: nêu cơ cấu 45 thật + 45 tổng hợp.
6. 5.1a: không gọi là "dưới lưới".
