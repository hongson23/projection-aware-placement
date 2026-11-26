# Bảng giờ chạy — lượt đo giờ ĐỘC QUYỀN (Bài 1 ch3, bản 3)

Bản v2 (25/11): mọi lượt DOC_QUYEN vào trung vị. So bit: lệch chỉ ở sai số làm tròn (≤ 1e-12): ['k08_CX2b_6512', 'k08_CX2e_6511', 'k08_CXe_6511', 'k08_CXe_6512', 'k08_CXe_6513', 'k08_P4_6511', 'k08_P4_6512', 'k08_P4_6513', 'k08_P8_6511', 'k08_P8_6512', 'k08_P8_6513', 'k08_U5b_6511', 'k08_U5b_6512', 'k08_U5b_6513', 'tv_real_nn1x_7311', 'tv_real_nn1x_7311_l2', 'tv_real_td_7311', 'tv_real_u5ktt_7311', 'tv_real_un1x_7311']; PHÂN KỲ (lệch lớn): [('k08_CX2b_6511', 0.0368), ('k08_CX2b_6513', 0.002), ('k08_CX2e_6512', 1.0), ('k08_CX2e_6513', 0.1881)].
(Dòng gốc bản 3:) Chỉ lượt DOC_QUYEN và không KHÁC bit. KHÁC bit (loại, phải giải thích): ['k08_P8_6511', 'k08_P4_6511', 'k08_U5b_6511', 'k08_P8_6512', 'k08_P4_6512', 'k08_U5b_6512', 'k08_P8_6513', 'k08_P4_6513', 'k08_U5b_6513', 'k08_CX2b_6511', 'k08_CX2e_6511', 'k08_CXe_6511', 'k08_CX2b_6512', 'k08_CX2e_6512', 'k08_CXe_6512', 'k08_CX2b_6513', 'k08_CX2e_6513', 'k08_CXe_6513', 'tv_real_td_7311', 'tv_real_nn1x_7311', 'tv_real_un1x_7311', 'tv_real_u5ktt_7311']. Thiếu/lỗi: []; không độc quyền (giữ, loại khỏi bảng): [].

## k08 (D5 khe 0,8), 3 seed × 1 lượt, hai làn song song — trung vị [min–max] trên seed

| nhánh | n | khởi động | khởi tạo | đặt | chuẩn bị PROJECT | PROJECT | mesh + đo | tổng phương pháp | HD (đánh giá) | tường | đỉnh GPU torch (GB) | đỉnh GPU nvidia-smi (GB) | đỉnh RSS (GB) | pp / P8 (theo seed) | so bit |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| U5b | 3 | 1.9 [1.9–1.9] | 13.8 [13.8–13.9] | 510.0 [473.7–519.5] | 5.4 [5.4–5.5] | 154.1 [136.1–181.7] | 8.9 [8.9–8.9] | 692.2 [637.9–729.5] | 0.4 [0.4–0.5] | 694.5 [640.3–731.8] | 16.03 [16.02–16.03] | 17.49 [17.47–17.55] | 2.17 [2.11–2.20] | 0.540 [0.503–0.559] | 0/3 TRUNG |
| P4 | 3 | 2.0 [1.9–2.1] | 13.8 [13.7–14.0] | 446.3 [405.6–463.1] | 5.4 [5.3–5.6] | 395.1 [383.0–427.9] | 9.1 [9.0–9.2] | 869.8 [862.0–874.3] | 0.4 [0.4–0.4] | 872.1 [864.5–876.7] | 16.04 [16.02–16.04] | 17.61 [17.53–17.65] | 2.19 [2.03–2.19] | 0.680 [0.648–0.702] | 0/3 TRUNG |
| P8 | 3 | 2.0 [1.9–2.0] | 14.2 [13.8–14.2] | 899.9 [820.7–936.3] | 5.6 [5.4–5.7] | 385.0 [311.0–419.1] | 9.0 [8.9–9.1] | 1268.4 [1239.2–1350.2] | 0.4 [0.4–0.4] | 1270.8 [1241.5–1352.6] | 16.04 [16.02–16.04] | 17.53 [17.51–17.65] | 2.07 [2.05–2.12] | 1.000 [1.000–1.000] | 0/3 TRUNG |
| CXe | 3 | 2.0 [2.0–2.1] | 13.9 [13.7–13.9] | 470.6 [428.6–480.2] | 5.2 [5.2–5.4] | 417.3 [406.3–429.8] | 9.2 [9.1–9.2] | 905.2 [886.5–925.9] | 0.4 [0.4–0.4] | 907.6 [889.0–928.3] | 16.04 [16.02–16.04] | 17.55 [17.40–17.61] | 2.11 [2.04–2.30] | 0.699 [0.686–0.730] | 0/3 TRUNG |
| CX2e | 3 | 2.0 [2.0–2.4] | 13.9 [13.7–14.0] | 504.6 [461.5–513.4] | 5.3 [5.3–5.3] | 439.5 [386.1–443.2] | 9.1 [9.1–9.2] | 932.8 [919.1–981.3] | 0.4 [0.4–0.4] | 935.2 [921.9–983.7] | 16.04 [16.02–16.04] | 17.51 [17.49–17.59] | 2.15 [2.13–2.19] | 0.735 [0.727–0.742] | 0/3 TRUNG |
| CX2b | 3 | 2.1 [1.9–2.2] | 14.0 [13.8–14.8] | 746.5 [728.5–762.1] | 5.1 [5.0–5.1] | 435.1 [356.0–488.5] | 9.1 [9.1–9.2] | 1225.1 [1131.5–1245.3] | 0.4 [0.4–0.4] | 1227.8 [1133.8–1247.8] | 16.04 [16.03–16.04] | 17.51 [17.51–17.53] | 2.12 [2.08–2.16] | 0.913 [0.907–0.982] | 0/3 TRUNG |

Ghi chú: CX* chạy bộ dựng v122u5x, U5b/P4/P8 chạy v122u5e (cùng phép tính khi cờ tắt); ở lượt cũ U5b trên u5x chậm hơn ≈ 3,5 % (u5_x_U5bcal 661/773 s so với 637/747 s) — tỉ số CX/P8 mang độ lệch bộ dựng này. Tầng 1 (U5b/P4/P8) mỗi seed trọn một làn; CX seed 6513 ở làn khác P8 6513.

## Thực vật TV-real seed 7311, 1 lượt — giờ từng bước do bộ dựng tự ghi

| nhánh | độc quyền | khởi tạo | đặt | PROJECT | tổng pp | HD | tường | đỉnh GPU torch (GB) | đỉnh RSS (GB) | igl phụ dùng | so bit | các bước (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| nn1x | DOC_QUYEN | 84.7 | 667.1 | 28.0 | 2396.0 | 266.3 | 2664.3 | 5.31 | 26.19 | [True] | KHAC | bo_tri 751.9, hop_con_psi 0.7, hang 0.1, chieu_dykstra 48.4, mc_bo_phan_ran 15.7, thong_ke_va_D_igl 1314.6, khac 1.3, vung_do_thi_ong 5.0, vung_chon_o 0.6, vung_truong_gpu 29.5, vung_mc_cpu 0.0, vung_ghep_trimesh 33.6, hd_vung_truoc 255.7, luu_npz 36.9, do_hang_ran_igl 0.0, do_vung_igl 28.4, do_gan_huong 0.0 |
| td | DOC_QUYEN | 85.3 | 627.6 | 37.4 | 2499.9 | 272.9 | 2774.9 | 5.54 | 36.44 | [True] | KHAC | bo_tri 712.9, hop_con_psi 0.8, hang 0.1, chieu_dykstra 59.5, mc_bo_phan_ran 15.7, thong_ke_va_D_igl 1443.2, khac 1.4, vung_do_thi_ong 5.0, vung_chon_o 0.6, vung_truong_gpu 25.6, vung_mc_cpu 0.0, vung_ghep_trimesh 34.5, hd_vung_truoc 259.5, luu_npz 38.2, do_hang_ran_igl 0.0, do_vung_igl 30.7, do_gan_huong 0.0 |
| u5ktt | DOC_QUYEN | 85.4 | 653.9 | 24.5 | 2287.5 | 268.1 | 2557.6 | 5.38 | 38.24 | [True] | KHAC | bo_tri 739.4, hop_con_psi 0.7, hang 0.1, chieu_dykstra 46.8, mc_bo_phan_ran 16.4, thong_ke_va_D_igl 1206.6, khac 1.4, vung_do_thi_ong 5.0, vung_chon_o 0.6, vung_truong_gpu 25.8, vung_mc_cpu 0.0, vung_ghep_trimesh 36.2, hd_vung_truoc 261.1, luu_npz 39.9, do_hang_ran_igl 0.0, do_vung_igl 31.9, do_gan_huong 0.0 |
| un1x | DOC_QUYEN | 83.9 | 524.6 | 26.8 | 2222.3 | 254.7 | 2478.9 | 5.35 | 35.20 | [True] | KHAC | bo_tri 608.6, hop_con_psi 0.7, hang 0.1, chieu_dykstra 47.5, mc_bo_phan_ran 15.3, thong_ke_va_D_igl 1289.6, khac 1.3, vung_do_thi_ong 4.5, vung_chon_o 0.6, vung_truong_gpu 22.1, vung_mc_cpu 0.0, vung_ghep_trimesh 33.0, hd_vung_truoc 244.3, luu_npz 37.2, do_hang_ran_igl 0.0, do_vung_igl 27.9, do_gan_huong 0.0 |

Cột igl phụ: cờ `dung` của tiến trình igl phụ; phải True ở mọi nhánh để giờ so được với nhau (lượt chính thực vật có khi chạy KHÔNG có igl phụ ⇒ bước D_igl chậm ≈ 2–3×; lượt đo luôn dùng igl phụ).
Tổng pp của thực vật = giay_tong − Σ hd_*: CÒN GỒM bước `thong_ke_va_D_igl` = D tới mesh sau chiếu (phương pháp: dùng cho vùng phụ thuộc) + D0 tới mesh TRƯỚC chiếu (chỉ đánh giá, khi có --hd_vung_truoc) + thể tích/IoU (đánh giá) — không tách được nếu không sửa mã; báo bước này riêng, và tổng pp là CẬN TRÊN của giờ phương pháp. Lượt --tv_thuan bỏ phần đánh giá.

## Nhiễu chéo hai làn

- Không đo (không bật --don hoặc thiếu lượt độc quyền): hai việc trên hai GPU dùng chung L3/bộ nhớ của máy; ghi là giới hạn.

## Độc quyền (giam_doc_quyen.py)

- 22/22 lượt DOC_QUYEN; chạy lại (_l2): 1; GPU lạ max 0; tính toán người dùng khác max 0.0000 lõi.
- Tải hệ thống (falco, tailscaled, kworker…; được phép) trên vùng lõi: trung vị 0.063 [0.027–0.081] lõi; lượt vượt nền p95: [].
- Tần số lõi việc TB 4006 [3788–4131] MHz; nhiệt GPU max 69.0 °C.

## So từng bit với lượt gốc

0/22 TRUNG; KHÁC: ['k08_P8_6511', 'k08_P4_6511', 'k08_U5b_6511', 'k08_P8_6512', 'k08_P4_6512', 'k08_U5b_6512', 'k08_P8_6513', 'k08_P4_6513', 'k08_U5b_6513', 'k08_CX2b_6511', 'k08_CX2e_6511', 'k08_CXe_6511', 'k08_CX2b_6512', 'k08_CX2e_6512', 'k08_CXe_6512', 'k08_CX2b_6513', 'k08_CX2e_6513', 'k08_CXe_6513', 'tv_real_td_7311', 'tv_real_nn1x_7311', 'tv_real_un1x_7311', 'tv_real_u5ktt_7311']; không có gốc: []. Thực vật không có gốc để so: [].

## Môi trường đo

```
CPU(s):                             32
Thread(s) per core:                 2
Model name:                         AMD Ryzen Threadripper PRO 5955WX 16-Cores
NUMA node0 CPU(s):                  0-31
index, name, driver_version, memory.total [MiB], pcie.link.gen.current, pcie.link.width.current
0, NVIDIA GeForce RTX 4090, 570.133.07, 24564 MiB, 1, 16
1, NVIDIA GeForce RTX 4090, 570.133.07, 24564 MiB, 1, 16
| NVIDIA-SMI 570.133.07             Driver Version: 570.133.07     CUDA Version: 12.8     |
python 3.11.16 torch 2.11.0+cu128 cuda 12.8 cudnn 91900
ondemand
```
