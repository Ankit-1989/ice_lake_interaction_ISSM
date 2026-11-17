# Proglacial lake filling and drainage: ISSM experiments

This repository holds the model code for the following paper:

> Pramanik, A., Greenwood, S. L., Morlighem, M., Barnett, J., Holmes, F. A., Gyllencreutz, R. and Regnell, C.:
> *Modelling the effects of proglacial lake filling and drainage on lake-terminating glacier dynamics.*
> Journal, volume, pages, DOI: `xx.xxxx/xxxxx` (to be added after publication)

The model is built with the [Ice-sheet and Sea-level System Model (ISSM)](https://issm.jpl.nasa.gov). It covers Qooqqup Sermia, a lake-terminating outlet glacier near Narsarsuaq in south Greenland. One MATLAB script, `runme_ice_lake.m`, builds the model, runs the inversion and spin-ups, and runs every experiment in the paper.

---

## Model summary

| Component        | Setting                                                                     |
|------------------|-----------------------------------------------------------------------------|
| Stress balance   | Shelfy-Stream Approximation (SSA)                                           |
| Basal sliding    | Budd law (p = q = 1), coefficient from inversion of observed surface velocity |
| Ice rheology     | Glen's flow law, B from Cuffey & Paterson (2010) at -5 °C                   |
| Surface mass balance | Constant field: MAR v3.11.3 monthly mean, 1960–2000                     |
| Grounding line   | Sub-element migration, no melt on partially floating elements               |
| Ice front        | Moving front with von Mises calving during spin-ups; fixed front in experiments |
| Lake             | Prescribed lake level on the lake nodes only; flat lake floor at -200 m     |
| Time step        | 1 month (1 hour for Exp 3b); output saved yearly (daily for Exp 3b)         |

### How the lake is represented

The lake level is applied through ISSM's sea-level field (`md.initialization.sealevel`), but only on the mesh nodes inside `Exp/lake_front_contour_points.exp`. Sea level is 0 m everywhere else. The lake floor is flattened to -200 m, and lake levels are given as **water depth above the lake floor**, as in the paper:

```
sealevel (lake nodes) = depth + lake_floor        e.g. a 100 m lake -> sealevel = -100 m
```

A lake level that changes over time is passed to ISSM as a forcing matrix. The last row holds the times, and ISSM interpolates linearly between them. Lake water density is 1000 kg m⁻³.

---

## Repository contents

```
.
├── runme_ice_lake.m        Main script: base model, spin-ups and all experiments
├── ice_lake.par            ISSM parameter file (geometry, velocities, rheology, ...)
├── Exp/
│   ├── NSQ_mouginot_modified_basin_w_LdFjLk.exp   Model domain outline
│   └── lake_front_contour_points.exp             Lake area where the lake level is applied
├── Scripts/                Helper functions (see Requirements)
├── Model_Data/             Input datasets (not included; see Input data)
└── Outputs/                Created by the script; model results are written here
```

---

## Requirements

- **ISSM**: built with MATLAB support. Installation guide: <https://issm.jpl.nasa.gov/download/>. Source: <https://github.com/ISSMteam/ISSM>.
- **MATLAB R2016b or newer**, because the script uses local functions.
- **Helper functions** in `Scripts/`, which the script loads with `addpath('./Scripts')`:
  - `interpBedmachineGreenland.m`: interpolates BedMachine Greenland onto the mesh
  - `interpMAR_monthly.m`: interpolates monthly MAR SMB onto the mesh
  - `interpSeaRISE_new.m`: interpolates SeaRISE fields (e.g. surface temperature) onto the mesh

### Input data

These datasets are not redistributed here. Download them from the original providers and put them in `Model_Data/` under the file names below.

| File in `Model_Data/`                 | Dataset                                                                 |
|---------------------------------------|-------------------------------------------------------------------------|
| `vx.nc`, `vy.nc`                      | MEaSUREs Greenland ice velocity mosaic (NSIDC), x and y components      |
| `BedMachineGreenland-v5.nc`           | BedMachine Greenland v5 (Morlighem et al.), <https://nsidc.org/data/idbmg4/versions/5> |
| `MARv3.11.3-historical-combined.nc`   | MAR v3.11.3 regional climate model, monthly surface mass balance        |
| `lake_mask.mat`                       | Lake outline mask used for the mesh refinement (included / provided by the authors) |

---

## Before running

1. Start MATLAB with ISSM on the path and `cd` into this folder.
2. In `ice_lake.par`, change the hard-coded `addpath` line to point to your copy of `Scripts/`.
3. Check the **User settings** block at the top of `runme_ice_lake.m`:
   - `cfg.np` is the number of MPI processes.
   - `cfg.smb_months` selects the MAR records that are averaged. The default `1:480` gives 1960–2000 if your MAR file starts in January 1960.
   - `nyrs_spinup` (20 yr) and `nyrs_exp` (100 yr) set the run lengths.

---

## Usage

The script uses ISSM's `organizer`. Choose the steps to run with the `steps` variable:

```matlab
steps = 1:4;    runme_ice_lake     % build the mesh and model, invert for basal friction
steps = 5:7;    runme_ice_lake     % spin-ups
steps = 8:15;   runme_ice_lake     % experiments
steps = 12;     runme_ice_lake     % any single step
```

Each step loads the saved output of the step it depends on from `Outputs/`, so steps can be run separately. The order must still be respected: run 1–4 before the spin-ups, and the spin-ups before the experiments.

### Steps

| Step | Name              | Description                                                               | Starts from |
|-----:|-------------------|---------------------------------------------------------------------------|-------------|
| 1    | `Mesh`            | Anisotropic mesh, refined on fast flow, the fjord and the lake            | –           |
| 2    | `Parameterization`| Geometry, masks, velocities and rheology (`ice_lake.par`)                 | 1           |
| 3    | `Geometry_update` | SSA, ice mask, flat lake floor at -200 m, boundary conditions             | 2           |
| 4    | `Inversion_Budd`  | Inversion for the Budd friction coefficient (M1QN3)                       | 3           |
| 5    | `SpinUp_Lake100`  | 20-yr spin-up with a 100 m lake (moving front)                            | 4           |
| 6    | `SpinUp_NoLake`   | 20-yr spin-up without a lake (moving front)                               | 4           |
| 7    | `SpinUp_Lake200`  | 20-yr spin-up with a 200 m lake, floating terminus (moving front)         | 4           |
| 8    | `Exp1_Lake100`    | **Exp 1**: 100 m lake held for 100 yr                                     | 5           |
| 9    | `Exp1_NoLake`     | **Exp 1**: land-terminating, 100 yr                                       | 6           |
| 10   | `Exp2a_Rise`      | **Exp 2a**: lake rises 100 → 200 m over 50 yr                             | 5           |
| 11   | `Exp2a_Fall`      | **Exp 2a**: lake falls 100 → 0 m over 50 yr                               | 5           |
| 12   | `Exp2b_Rise`      | **Exp 2b**: lake rises 0 → 100/125/150/175/200 m over 50 yr               | 6           |
| 13   | `Exp2b_Fall`      | **Exp 2b**: lake falls 100/125/150/175/200 → 0 m over 50 yr               | 12 (year 70) |
| 14   | `Exp3a_Drainage`  | **Exp 3a**: 100 m drop over 1, 2, 3 or 4 yr, from a 200 m and a 100 m lake | 7 and 5    |
| 15   | `Exp3b_GLOF`      | **Exp 3b**: 10/20/30/40/50 m drop in 1 day, from a 200 m lake (hourly step) | 7         |

Steps 8–15 all use a fixed ice front.

**Run times.** The base model and the yearly experiments take hours. Exp 3b uses an hourly time step and is much more expensive per model year. To check the whole chain quickly, first run with `nyrs_spinup = 1` and `nyrs_exp = 2`.

---

## Outputs

All results are written to `Outputs/` as MATLAB `.mat` files:

- Single-run steps are saved by the organizer as `Outputs/Narsarsuaq<StepName>.mat`, e.g. `NarsarsuaqExp1_Lake100.mat`. Load one with `md = loadmodel(org,'Exp1_Lake100')`.
- Steps with several runs save one file per run with a descriptive name:

| Step | File name pattern                                   | Example                                  |
|-----:|-----------------------------------------------------|------------------------------------------|
| 12   | `Narsarsuaq_Exp2b_Rise_0to<D>.mat`                  | `Narsarsuaq_Exp2b_Rise_0to150.mat`       |
| 13   | `Narsarsuaq_Exp2b_Fall_<D>to0.mat`                  | `Narsarsuaq_Exp2b_Fall_150to0.mat`       |
| 14   | `Narsarsuaq_Exp3a_from<D>_drop100_in<N>yr.mat`      | `Narsarsuaq_Exp3a_from200_drop100_in2yr.mat` |
| 15   | `Narsarsuaq_Exp3b_GLOF_drop<dd>m_1day.mat`          | `Narsarsuaq_Exp3b_GLOF_drop30m_1day.mat` |

Results are in `md.results.TransientSolution`. Each entry includes ice thickness, velocity, the grounded and floating areas, grounding-line flux, the ocean and ice level-sets, the von Mises stress and the deviatoric stresses.

---

## Citation

If you use this code, please cite the paper above and this archive:

> Pramanik, A. et al. (2026). Code for "Modelling the effects of proglacial lake filling and drainage on lake-terminating glacier dynamics" (v1.0). Zenodo. DOI: `10.5281/zenodo.xxxxxxx`

Please also cite ISSM (Larour et al., 2012, *J. Geophys. Res.*, doi:10.1029/2011JF002140) and the input datasets listed above.

## License

To be added (e.g. MIT, BSD-3-Clause or CC-BY-4.0). ISSM itself is distributed under its own BSD-style license.

## Contact

Ankit Pramanik: please open an issue in this repository for questions about the code.
