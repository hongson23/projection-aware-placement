# Environment used for every number in the paper

All builder, checker and analysis runs used one conda environment (internally called `b2demo`) on one machine.

| item | value |
|---|---|
| OS | Ubuntu 20.04.5 LTS, Linux 5.4.0 |
| CPU | AMD Ryzen Threadripper PRO 5955WX (16 cores) |
| GPU | 2 x NVIDIA GeForce RTX 4090 (24 GB each), driver 570.133.07, CUDA 12.8 (driver) |
| Python | 3.11.16 |
| PyTorch | 2.11.0+cu128 |
| key packages | numpy 2.4.6, scipy 1.17.1, trimesh 5.1.0, libigl 2.6.3, scikit-image 0.26.0, manifold3d 3.5.4, tifffile 2026.3.3 |
| full list | `requirements_b2demo.txt` (`python -m pip freeze`, captured read-only on 25/11/2026) |

The environment also contains packages that the code in this repository does not import (for example SAM-2 and diffusers); they come from other projects that share the environment. A minimal install for the builder and checker is

```
pip install torch==2.11.0 numpy==2.4.6 scipy==1.17.1 trimesh==5.1.0 libigl==2.6.3 scikit-image==0.26.0 manifold3d==3.5.4 tifffile==2026.3.3
```

but only the frozen list reproduces the runs exactly.

## External baselines

* cellPACK ran in a separate environment (Python 3.11.18) with cellPACK installed in editable mode from
  `mesoscope/cellpack` at commit `a576c5c1e81a16d7b908bfba9168e708f4730ce7`; its frozen package list is
  `requirements_cellpack.txt`. The GPU port is a patch on top of that commit (`external_baselines/cellpack/gpu_patch/`);
  it imports torch from the main environment through the variable `CELLPACK_TORCH_SITE`.
* Packmol: official release 100.0.0 (`packmol-100.0.0-linux-x86_64.tar.gz`); the release checksums are in
  `external_baselines/packmol/SHA256SUMS_packmol-100.0.0.txt`. The binary is not redistributed here.

## Rendering

The 3D figures were rendered with Blender 4.2.1 LTS (Cycles, GPU); see `figures/HINH_CAD.md`.
