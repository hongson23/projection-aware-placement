# Đối chiếu luật sự thật v1.2 với bộ kiểm độc lập (11/11/2026)

Trạng thái: **ĐÓNG BĂNG 11/11/2026 — đã qua 2 vòng rà soát chéo độc lập**. Sinh tự động bởi `tao_luat_su_that.py` (cùng lượt với ba tệp JSON).

Nguồn số liệu bộ kiểm (chỉ đọc, không sửa gì trên server): `~/Vision_team_data_189/sonhh20/ch3_tebao_dep/bo_kiem_v2/BAO_CAO.md` (động vật, 05/11, đã qua 2 vòng rà nhưng bản sau M.h chưa rà lại) và `BAO_CAO_THUCVAT.md` (thực vật, 06/11, **CHƯA rà soát chéo**; người viết tự ghi là đã đọc Phụ lục 2c nên không độc lập hoàn toàn với trường máy). Mô tả CHECK: `Inside3D/TaiLieu_PhaiNho.md` mục 7.

Cột:
- **Bộ kiểm**: có / có nhưng kiểm giá trị tham số (CO*) / một phần / không kiểm được trên mesh hiện có / chưa có mục riêng.
- **Đột biến**: số đột biến nhắm luật đó bị chính luật đó bắt / tổng (bảng độ phủ của báo cáo).
- **Tin cậy**: trong danh sách tin cậy anh Sơn đưa (chỉ có cho động vật). Thực vật: chưa có danh sách vì bộ kiểm chưa rà.
- **Phụ thuộc ngưỡng tự đặt**: id trong `tham_so_quy_uoc_v1_2_20261111.json`.
- **Nhóm A/B/C** (theo rà soát chéo 11/11): **A** = bộ kiểm đo giá trị vẽ/chọn (hoặc dải cũ) thay cho dải sự thật; **B** = ngưỡng của bộ kiểm chính là một tham số mà luật sự thật cần để định lượng; **C** = v1 xếp nhầm.
- **Lưu ý thực vật:** bộ kiểm thực vật mới chỉ thử trên 3 tế bào do chính bộ kiểm dựng (TV1, TV2, TV6), chưa có mesh từ bộ dựng thật.
- Phán quyết trong bảng là của báo cáo bộ kiểm theo luật v4.1/v1.4; **bộ kiểm CHƯA chạy lại trên v1.2** (các luật đổi dải/vị từ: 4.3b, 1.2a, 2.1, 2.2a, 2.4b, 10.2b, 10.2c, 5.4a, 5.6a, 10.3b, 4.5b, 4.6b, 4.8b).

## Tế bào động vật (20 luật sự thật)

| Luật | Lời (sự thật) | Bộ kiểm | Cách kiểm hiện có | Nhóm A/B/C | Đột biến | Tin cậy | Phụ thuộc ngưỡng tự đặt |
|---|---|---|---|---|---|---|---|
| 1.1 | Tế bào có ba thành phần: màng tế bào bao ngoài, chất tế bào bên trong, nhân. | có | đúng 1 màng + 1 nhân, chất tế bào > 0, nhân trong màng | — | 2/2 | có | `chung.eps_giao_um3` |
| 1.2a | Tế bào là khối ba chiều; tế bào người tròn/trôi nổi đường kính cỡ 10 µm trở lên (HeLa trôi nổi 17,1 µm; thể tích HeLa khoảng 500–4 000 µm³). | có | 3 cạnh OBB của màng ∈ [10; 30] µm — cận trên 30 là tham số chọn | B: cận trên 30 µm = dv.can_tren_duong_kinh_te_bao | 2/2 | có (phần cỡ) | `dv.can_tren_duong_kinh_te_bao` |
| 1.3 | Thường có một nhân, nằm trong khối chất tế bào, không bắt buộc ở giữa. | có | 1 nhân, liền khối, trong tế bào; lệch tâm chỉ ghi thông tin | — | 2/2 | có | — |
| 1.4a | Nhân nằm trong chất tế bào. | có | bool nhân − màng rỗng | — | 2/2 | có | `chung.eps_giao_um3` |
| 1.5 | Đường kính nhân bằng khoảng 0,35–0,6 đường kính tế bào; nhân rộng 5–10 µm. | có | (i) tỉ số đường kính cầu cùng thể tích ∈ [0,35; 0,6]; (ii) mọi cạnh OBB nhân ∈ [5; 10] | — | 2/2 | có (phần cỡ) | `chung.dung_sai_khoang` |
| 4.2a | Trong nhân có 1–5 nhân con, mỗi cái rộng 0,5–5 µm, không có màng bao. | một phần | kiểm SỐ VẼ 1–3 (tham số dv.so_nhan_con) chứ không 1–5; cỡ [0,5; 5]; trong nhân. 'Không có màng bao' không kiểm được | A: đếm 1–3 (vẽ) thay 1–5 | 6/6 | có (phần cỡ) | — |
| 4.3a | Ti thể nằm trong chất tế bào. | có | trong màng và không giao khối đặc nhân | — | 4/4 | có | `chung.eps_giao_um3` |
| 4.3b | Ti thể có đường kính tiết diện (độ dày ống/hạt) 0,2–1 µm. | có | rộng = độ dày cục bộ đo theo dải CŨ [0,5; 1] µm (v4.1), chặt hơn sự thật v1.1 [0,2; 1]; dài [1; 2] nay là tham số dv.dai_ti_the | A: rộng [0,5; 1] (cũ) thay [0,2; 1] | 10/10 | có | `kiem.dv.NG_CONG_43B`, `dv.dai_ti_the` |
| 4.4a | Có lưới nội chất hạt (có ribosome) và lưới nội chất trơn; lưới nội chất liền với màng ngoài của nhân và trải khắp chất tế bào (dạng tấm quanh nhân, dạng ống ở xa). | một phần | kiểm ER chạm mặt ngoài nhân (≤ 0,01 µm), khớp vị từ 'chạm, diện tích > 0'; không kiểm ER hạt/trơn, 'trải khắp' | B: 0,01 µm = chung.tiep_xuc_toi_da_um | 2/2 | không | `chung.tiep_xuc_toi_da_um` |
| 4.5a | Bộ máy Golgi gồm 4–11 túi dẹt xếp chồng, mỗi túi rộng khoảng 0,7–1,1 µm. | không kiểm được trên mesh hiện có | Golgi là MỘT khối ở cả 6 tế bào → không đếm/đo được túi; nếu tách được thì kiểm số vẽ 4–8 | — | không có đột biến | không | — |
| 4.5b | Golgi nằm cạnh nhân, ôm quanh trung thể. | có | định nghĩa tương đối của bộ kiểm = vị từ sự thật v1.1 (khe tới nhân < khe tới màng; tâm trung thể trong bao lồi Golgi); luật mềm | — | không có đột biến | không | — |
| 4.6a | Lysosome là túi nhỏ đường kính 0,1–1,2 µm (thường 0,1–0,5 µm). | có | mọi cạnh OBB ∈ [0,1; 1,2] | — | 2/2 | có (phần cỡ) | — |
| 4.6b | Lysosome tập trung quanh nhân và rải ra phía ngoài. | có | trung vị u < 0,5 và có u ≥ 0,5 = vị từ sự thật v1.1; luật mềm | — | không có đột biến | không | — |
| 4.8a | Trung thể gồm hai trung tử đặt vuông góc, mỗi trung tử dài 0,3–0,5 µm, rộng 0,15–0,25 µm. | không kiểm được trên mesh hiện có | trung thể là MỘT khối ở 5/6 tế bào → không đo được trung tử; ≥ 3 khối = TRƯỢT | — | 2/2 | không | `dv.goc_vuong_toi_thieu` |
| 4.8b | Trung thể nằm gần nhân. | có | khe trung thể–nhân < khe trung thể–màng = vị từ sự thật v1.1; luật mềm | — | 2/2 | không | — |
| M.a | Mọi bào quan nằm trong màng tế bào. | có | mọi bộ phận − khối màng ≤ ε; bộ phận mất = trượt | — | 4/4 | có | `chung.eps_giao_um3` |
| M.d | Màng nhân gồm hai lớp, khoang giữa 20–50 nm, có lỗ nhân. | không kiểm được trên mesh hiện có | chỉ kiểm được khi nhân là vỏ rỗng; khi đó kiểm khoang VẼ (≤ 0,3 µm); số thật 20–50 nm không kiểm (đúng ý vị từ) | — | 2/2 | không | `kiem.dv.MD_KHOANG_MAX`, `kiem.dv.MD_TI_LE_TIA` |
| M.e | Màng trong ti thể gấp thành mào. | không kiểm được trên mesh hiện có | mesh ti thể là mặt bao; không có màng trong/mào | — | không có đột biến | không | — |
| M.c₀ | Hai bào quan khác loại không giao nhau (khe ≥ 0; được chạm). | chưa có mục riêng | không có mục riêng. Phép bool giao của M.c (khe 0,3) bao luôn 'không giao' cho mọi cặp KHÔNG miễn; cặp miễn MIEN_MC (ER–nhân, ribosome–ER) hiện không được kiểm giao | — | — | suy từ M.c (tin cậy), chưa tách | `chung.eps_giao_um3`, `kiem.dv.MIEN_MC` |
| M.h | Hai bộ phận cùng loại không chồng lên nhau (được chạm, không được giao). | có | bool giao cùng loại > 1e-6 µm³ = trượt; chạm được phép | — | 3/3 | có | `chung.eps_giao_um3` |

## Tế bào thực vật (28 luật sự thật)

| Luật | Lời (sự thật) | Bộ kiểm | Cách kiểm hiện có | Nhóm A/B/C | Đột biến | Tin cậy | Phụ thuộc ngưỡng tự đặt |
|---|---|---|---|---|---|---|---|
| 2.1 | Tế bào thực vật có thành tế bào ở ngoài cùng; bên trong là màng tế bào, chất tế bào, nhân; trong chất tế bào có không bào trung tâm và lục lạp. | có | có đủ thành phần, nằm đúng chỗ; bộ kiểm đòi đúng 1 nhân, vị từ v1.1 chỉ đòi ≥ 1 | — | 1/1 | chưa có danh sách | `chung.eps_giao_um3` |
| 2.2a | Tế bào thịt lá trưởng thành có thể tích phía trong màng tế bào (phần sống) 60–100 pL (mô xốp ~60 pL, mô giậu ~100 pL). | có, nhưng kiểm giá trị THAM SỐ (vẽ/tự đặt) | kiểm GIÁ TRỊ CHỌN cũ (hộp ngoài 60×30×30 ± 2 %, 54 pL ± 10 %) — KHÔNG kiểm thể tích trong màng 60–100 pL; hộp v1.2 đã đổi | A: hộp ngoài 60×30×30 (chọn, cũ) thay thể tích trong màng 60–100 pL | 1/1 | chưa có danh sách | `kiem.tv.TV_TOL_HOP`, `kiem.tv.TV_TOL_THE_TICH_HOP`, `tv.kich_thuoc_vach` |
| 2.3a | Thành tế bào nằm ngoài màng tế bào và bao quanh toàn bộ tế bào. | có | màng trong khối thành, không giao vật liệu thành | — | 1/1 | chưa có danh sách | `chung.eps_giao_um3` |
| 2.4a | Tế bào trưởng thành có một không bào trung tâm, có màng riêng bao quanh. | một phần | đúng 1 không bào, 1 khối; 'màng riêng' không kiểm được | — | 1/1 | chưa có danh sách | — |
| 2.4b | Không bào trung tâm chiếm 68–85 % thể tích phía trong màng tế bào. | có | V(không bào)/V(trong màng) ∈ [0,70; 0,85] (v1.4) — chặt hơn sự thật v1.1 [0,68; 0,85] | A: [0,70; 0,85] (v1.4) thay [0,68; 0,85] | 1/1 | chưa có danh sách | — |
| 2.6b | Nhân nằm trong lớp chất tế bào sát thành (thường ở mặt đáy tế bào; dời ra thành bên khi chiếu sáng mạnh), không ở giữa không bào. | có, nhưng kiểm giá trị THAM SỐ (vẽ/tự đặt) | khe nhân–màng ∈ [0,3; 0,5] (cận dưới 0,3 là khe quy ước) và không giao không bào | A: cận dưới 0,3 µm là khe quy ước | 1/1 | chưa có danh sách | `tv.khe_nhan_mang`, `tv.gan_nhan_mang_toi_da` |
| 2.7b | Lục lạp nằm trong lớp chất tế bào sát màng tế bào, nằm rải trên mặt trong tế bào (không trôi giữa không bào); hướng mặt dẹt ở 2.7f. | có, nhưng kiểm giá trị THAM SỐ (vẽ/tự đặt) | khe tới màng ∈ [0,3; 0,5] (cận dưới 0,3 là khe quy ước); 'rải' ≥ 3 mặt hộp = vị từ 'rai' v1.1 | A: cận dưới 0,3 µm là khe quy ước | 3/3 | chưa có danh sách | `tv.khe_luc_lap_mang`, `tv.so_mat_rai` |
| 2.7c | Lục lạp hình thấu kính (đĩa dẹt hai mặt lồi). | có | khối lồi, đĩa tròn, hồ sơ thấu kính | — | 1/1 | chưa có danh sách | `kiem.tv.TV_NG_LOI`, `kiem.tv.TV_NEP_SAU`, `kiem.tv.TV_NG_TRON_DIA`, `kiem.tv.TV_THAU_KINH` |
| 2.7d | Lục lạp thực vật bậc cao đường kính khoảng 4–10 µm, dày khoảng 1,1–2,5 µm. | có, nhưng kiểm giá trị THAM SỐ (vẽ/tự đặt) | kiểm dải VẼ 5–7 × 2–2,5 (tham số) ⊂ dải sự thật 4–10 × 1,11–2,5 ⇒ đạt bộ kiểm thì đạt sự thật, nhưng trượt bộ kiểm chưa chắc trái sự thật | A: 5–7 × 2–2,5 (vẽ) thay 4–10 × 1,11–2,5 | 2/2 | chưa có danh sách | `tv.kich_thuoc_luc_lap_ve` |
| M.d | Ti thể, Golgi, peroxisome nằm trong lớp chất tế bào (sát màng tế bào hoặc sát nhân), không trôi giữa không bào. | có, nhưng kiểm giá trị THAM SỐ (vẽ/tự đặt) | cách màng hoặc nhân ≤ 1 µm (tham số); không trôi giữa không bào suy từ khe K.1 | B: 1 µm = tv.gan_mang_hoac_nhan | 2/2 | chưa có danh sách | `tv.gan_mang_hoac_nhan` |
| M.h | Hai bộ phận cùng loại không chồng lên nhau (được chạm, không được giao). | có | mọi cặp cùng loại (cả thylakoid trong granum, túi trong chồng): giao ≤ ε | — | 2/2 | chưa có danh sách | `chung.eps_giao_um3` |
| 5.1a | Trong nhân có một hoặc nhiều nhân con (hạch nhân), không có màng bao. | một phần | kiểm số VẼ 1–2 (tham số) và nằm trong nhân; 'không màng bao' không kiểm được | A: đếm 1–2 (vẽ) thay ≥ 1 | 1/1 | chưa có danh sách | `tv.so_nhan_con` |
| 5.2a | Màng nhân gồm hai lớp, giữa là khoang hẹp (thật 20–50 nm), có nhiều lỗ nhân; màng ngoài nối liền lưới nội chất. | có | ≥ 50 % tia gặp hai lớp, có khoang, có lỗ, chạm ER; khoang thật chỉ là thông tin (đúng ý vị từ) | C: v1 xếp nhầm vào 'kiểm tham số' (CO*); thực ra kiểm cấu trúc sự thật, khoang thật chỉ là thông tin → CO | 1/1 | chưa có danh sách | `kiem.tv.TV_MD_TIA` |
| 5.3e | Màng trong ti thể gấp nếp thành mào. | có, nhưng kiểm giá trị THAM SỐ (vẽ/tự đặt) | kiểm 4–8 mào (tham số) ⇒ có mào | A: đếm 4–8 mào (vẽ) thay ≥ 1 | 2/2 | chưa có danh sách | `tv.mao_ve` |
| 5.4a | Có lưới nội chất hạt (có ribosome) và trơn; màng lưới nội chất nối liền màng ngoài của nhân. | một phần | diện tích chạm ER–màng nhân ≥ 1 µm² (tham số); không phân biệt ER hạt/trơn | A: ≥ 1 µm² (tham số tv.vung_noi_er_mang_nhan) thay 'diện tích > 0' | 1/1 | chưa có danh sách | `tv.vung_noi_er_mang_nhan`, `chung.tiep_xuc_toi_da_um` |
| 5.6a | Bộ máy Golgi ở tế bào thực vật gồm nhiều thể Golgi (chồng túi dẹt) rời nhau, rải khắp chất tế bào; không gom thành một dải cạnh nhân. | có, nhưng kiểm giá trị THAM SỐ (vẽ/tự đặt) | chồng cách nhau ≥ 2 µm (tham số); có chồng cách nhân > 2 µm | A: 2 µm (tham số tv.khe_golgi_golgi) thay 'khe > 0' | 2/2 | chưa có danh sách | `tv.khe_golgi_golgi`, `kiem.tv.TV_GOLGI_GAN_NHAN` |
| 5.7a | Peroxisome là túi nhỏ một màng, dạng tròn (cầu hoặc bầu dục), đường kính 0,1–1 µm. | có, nhưng kiểm giá trị THAM SỐ (vẽ/tự đặt) | kiểm dải VẼ [0,5; 1] (tham số) ⊂ sự thật [0,1; 1] + độ tròn | A: [0,5; 1] (vẽ) thay [0,1; 1] | 1/1 | chưa có danh sách | `tv.duong_kinh_peroxisome_ve`, `kiem.tv.TV_NG_CAU` |
| 5.9 | Không có trung thể (trung tử). | có | không có bộ phận 'trung_the' | — | 1/1 | chưa có danh sách | — |
| 10.1a | Lục lạp có hai màng (màng ngoài, màng trong) bao quanh chất nền; chất nền có DNA và ribosome riêng. | một phần | một màng trong, khe 0,05 ± 0,01 (tham số); DNA/ribosome chất nền không vẽ | A: khe 0,05 ± 0,01 (vẽ) thay 'chứa, t = 0' | 1/1 | chưa có danh sách | `tv.khe_mang_trong_lap`, `kiem.tv.TV_TOL_KHE_LAP` |
| 10.2a | Trong chất nền lục lạp có các grana; mỗi granum là một chồng thylakoid dẹt. | có, nhưng kiểm giá trị THAM SỐ (vẽ/tự đặt) | kiểm 10–20 grana (tham số); chồng dẹt, song song, đồng trục | A: đếm 10–20 grana (vẽ) thay ≥ 1 | 1/1 | chưa có danh sách | `tv.so_granum`, `kiem.tv.TV_PHANG`, `kiem.tv.TV_GOC_THANG`, `kiem.tv.TV_LECH_CHONG` |
| 10.3a | Các grana nối với nhau bằng tấm thylakoid chất nền: toàn bộ hệ thylakoid là một khối màng liền. | có | đồ thị granum–tấm nối liền một thành phần trong mỗi lục lạp | — | 1/1 | chưa có danh sách | — |
| M.a | Mọi bào quan nằm trong màng tế bào. | có | mọi bộ phận trên cùng trong khối trong màng; con trong cha | — | 2/2 | chưa có danh sách | `chung.eps_giao_um3` |
| 2.6e | Nhân chiếm 0,16–2,3 % thể tích phía trong màng. | có | V(nhân)/V(trong màng) ∈ [0,16 %; 2,3 %] | — | 1/1 | chưa có danh sách | — |
| 10.2b | Granum đường kính khoảng 0,2–0,6 µm. | có | hai cạnh lớn OBB ∈ [0,3; 0,6] (v1.4) — chặt hơn sự thật v1.1 [0,2; 0,6] | A: [0,3; 0,6] (v1.4) thay [0,2; 0,6] | 1/1 | chưa có danh sách | — |
| 10.2c | Mỗi granum gồm từ 3 tới khoảng 100 thylakoid (thường 3–19, trung bình ~9 ở rau diếp thích nghi tối; dưới 10 ở cây ưa sáng mạnh, tới 100 ở cây rất ưa bóng). | có, nhưng kiểm giá trị THAM SỐ (vẽ/tự đặt) | kiểm số VẼ 3–8 (tham số) ⊂ sự thật v1.1 3–100 | A: đếm 3–8 (vẽ) thay 3–100 | 1/1 | chưa có danh sách | `tv.so_thylakoid_ve` |
| 2.7f | Lục lạp đặt mặt dẹt song song với màng tế bào. | có | trục ngắn PCA lệch pháp tuyến màng ≤ 30°; luật mềm | B: 30° = tv.lech_huong_luc_lap | 1/1 | chưa có danh sách | `tv.lech_huong_luc_lap` |
| 10.3b | Mỗi granum có ít nhất một tấm nối chạm vào nó, phải còn đúng sau bước chiếu. | có | mỗi granum có tấm nối chạm, diện tích ≥ 0,005 µm² (tham số) | A: ≥ 0,005 µm² (tham số tv.vung_noi_tam_granum) thay 'diện tích > 0' | 2/2 | chưa có danh sách | `tv.vung_noi_tam_granum` |
| M.c₀ | Hai bộ phận khác loại không giao nhau (khe ≥ 0; được chạm). | chưa có mục riêng | không có mục riêng. K.1 (khe 0,3/0,1) bao 'không giao' cho cặp KHÔNG miễn; 9 cặp miễn (ER–nhân, ER–màng nhân, ribosome–ER, tấm nối–granum, …) chưa được kiểm giao | — | — | chưa có danh sách | `chung.eps_giao_um3`, `tv.mien_khe_chung` |

## Luật quy ước cũ của động vật: nay là phép kiểm THAM SỐ (không phải luật)

| Luật cũ | Bộ kiểm | Đột biến | Tin cậy | Tham số |
|---|---|---|---|---|
| 1.2b | có | 2/2 | không | `kiem.dv.NG_12B_TRON` |
| 1.4b | có | 2/2 | có | `dv.T_nhan_mang` |
| 1.6 | KHÔNG ÁP DỤNG (Sinh10) | không có đột biến | không | `dv.ti_the_lop6` |
| 4.2b | có | 2/2 | có | `dv.T_nhancon` |
| 4.3c | có | 1/1 | có | `dv.so_ti_the` |
| 4.3d | có | 14/14 | có | `chung.khe_roi, chung.khe_cung_loai` |
| 4.4b | có | 2/2 | không | `dv.D_er, kiem.dv.DUNG_SAI_44B` |
| 4.5c | có | 2/2 | có | `chung.khe_roi` |
| 4.5d | MP (đã bỏ) | không có đột biến | không | `kiem.dv.NG_45D_TI_LE` |
| 4.6c | có | 1/1 | có | `dv.so_lysosome` |
| 4.6d | có | 2/2 | không | `kiem.dv.NG_46D_CAU` |
| 4.7 | không kiểm được trên mesh hiện có | không có đột biến | không | `dv.ribosome_ve` |
| 4.8c | có | 4/4 | có | `chung.khe_roi` |
| M.b | có | 2/2 | có | `dv.T_mang` |
| M.c | có | 12/12 | có | `chung.khe_roi` |
| M.f | có | 12/12 | có | `chung.lien_khoi_moi_bao_quan` |

## Tóm tắt

- Không có phép kiểm (không kiểm được trên mesh hiện có hoặc chưa có mục riêng): dv 4.5a, dv 4.8a, dv M.d, dv M.e, dv M.c₀, tv M.c₀.
- Chỉ kiểm một phần: dv 4.2a, dv 4.4a, tv 2.4a, tv 5.1a, tv 5.4a, tv 10.1a.
- Bộ kiểm đo GIÁ TRỊ THAM SỐ thay cho dải sự thật (cần tách phán quyết sự thật / tham số): dv 4.2a, tv 2.2a, tv 2.6b, tv 2.7b, tv 2.7d, tv M.d, tv 5.3e, tv 5.6a, tv 5.7a, tv 10.2a, tv 10.2c.
- Luật mềm động vật 4.5b, 4.6b, 4.8b: bộ kiểm dùng định nghĩa tương đối, bộ dựng dùng ngưỡng tuyệt đối (1,0 / 4,0 / 1,5 µm) — hai bên đang đo hai thứ khác nhau.
- Mọi kết luận 'đạt/trượt' trong bảng là số của báo cáo bộ kiểm, chưa được tính lại ở đây.

## Câu hỏi mở

1. 2.2a đã chốt (PI 11/11): thể tích TRONG màng; hộp và nhân v1.2 theo KKT_V12. Bộ kiểm 2.2a vẫn đo hộp ngoài cũ — cần sửa.
2. Tách vế sự thật từ luật nhãn quy ước (5.1a, 5.3e, 10.2a thực vật): PI chưa duyệt.
3. 4.4a: vế 'ER trải khắp chất tế bào' là sự thật nhưng mô hình không thể hiện (ghi ở khong_the_hien_co_chu_dich). Cách nói trong bài?
4. Bộ kiểm phải chạy lại với dải v1.1 (4.3b 0,2–1; 2.4b 0,68–0,85; 10.2b 0,2–0,6; 10.2c 3–100; 2.2a 60–100 pL; 2.1 ≥ 1 nhân) và tách phán quyết sự thật / tham số (nhóm A).
5. M.c₀: bộ kiểm chưa có mục riêng; cặp miễn khe (ER–nhân, ribosome–ER ở dv; 9 cặp K.1 ở tv) chưa được kiểm giao.
6. Mặc định botri2 còn trái sự thật (mau_thuan_mac_dinh_voi_su_that): trung thể một cầu (4.8a), Golgi một khối (4.5a), mien_cung_loai = lysosome (M.h), nhân đặc (M.d), ti thể không mào (M.e). Sửa bộ dựng hay ghi là giới hạn?
7. 10.2b cận dưới 0,2 µm ≈ trung bình − 2 SD của Bussi (tự tính); 10.2c cận trên 100 lấy từ câu tổng quan Staehelin (cây cực ưa bóng) — phạm vi có quá rộng cho tế bào lá thường không?
8. 2.6a đã chuyển thành tham số (không có nguồn khoa học nói thẳng 'một nhân' trong ghi chú) — cần tra thêm nguồn nếu muốn giữ là luật.
