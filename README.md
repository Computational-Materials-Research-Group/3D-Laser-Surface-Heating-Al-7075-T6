# 3D Laser Surface Heating — Al 7075-T6
### FreeFEM++ Adaptive Moving Mesh Simulation

![HAZ field showing concentric heat-affected rings around the laser spot](preview.png)

---

## Overview

This repository contains a **FreeFEM++** simulation of a high-power laser traversing a thin aluminium sheet (Al 7075-T6, 1 mm thick). The code solves the **3D transient heat conduction equation** with a moving Gaussian–Beer-Lambert laser source on a fully adaptive tetrahedral mesh that is rebuilt and refined around the laser spot at every time step.

The primary outputs are:

- Full 3D temperature field evolving over time
- Heat-Affected Zone (HAZ) boundary
- Melt pool extent (liquid aluminium region)
- Recast layer (re-solidified melt)
- Burn-through depth tracking per step
- Thermal gradient field (residual stress driver)

---
<img width="1008" height="772" alt="laser surf" src="https://github.com/user-attachments/assets/f9f62661-ef16-46eb-a1f4-1506589fb9ea" />

## Physics

### Governing equation

```
rho * cp * dT/dt  -  div( kc * grad(T) )  =  Q(x, y, z, t)
```

Discretised by **backward Euler** (implicit, unconditionally stable).

### Laser source model

```
Q(x,y,z,t) = Q0 * exp( -[(x-xc)² + (y-yc)²] / (2·σ²) )
                * exp( α·z )
```

| Term | Description |
|---|---|
| `Q0` | Peak volumetric power density [W/m³] |
| `σ` | Gaussian beam sigma (beam radius ≈ 2σ) |
| `α` | Beer-Lambert absorption depth coefficient [1/m] |
| `z ∈ [-Lz, 0]` | z=0 is top surface; source decays with depth |

### Boundary conditions

| Face | Condition |
|---|---|
| All four side faces | Dirichlet: T = T_amb (large plate assumption) |
| Top face (z = 0) | Laser source enters via volumetric term |
| Bottom face (z = −Lz) | Heat conducted through thickness |

---

## Material Properties — Al 7075-T6

| Property | Symbol | Value | Units |
|---|---|---|---|
| Thermal conductivity | kc | 130 | W/(m·K) |
| Density | ρ | 2810 | kg/m³ |
| Specific heat | cp | 960 | J/(kg·K) |
| Ambient temperature | T_amb | 298 | K (25 °C) |
| Solidus | T_sol | 750 | K (477 °C) |
| Liquidus | T_liq | 908 | K (635 °C) |
| Vaporisation | T_vap | 2792 | K (2519 °C) |
| HAZ threshold | T_HAZ | 573 | K (300 °C) |

---

## Adaptive Mesh Strategy

The mesh is **rebuilt every time step** — this is the key feature of the code.

```
Each step:
  1. Build coarse 2D background mesh (60×40 elements)
  2. adaptmesh × Nadapt passes driven by Gaussian indicator at (xc, yc)
       → fine elements concentrate at the laser spot
       → coarse elements far from the beam
  3. buildlayers → extrude adapted 2D mesh to 3D tetrahedra (nz layers)
  4. Interpolate T_old from fixed background grid onto new mesh
       → clamp any NaN to T_amb (robustness guard)
  5. Solve heat equation on adapted 3D mesh (UMFPACK direct solver)
  6. Save solution back to background grid for next step
```

### Mesh parameters

| Parameter | Value | Description |
|---|---|---|
| `Nadapt` | 3 | Adaptation passes per step |
| `err` | 0.004 | Interpolation error target |
| `hmin` | σ/20 = 0.025 mm | Finest element near beam centre |
| `hmax` | Lx/25 = 4 mm | Coarsest element far from beam |
| `nbvx` | 120 000 | Max vertices in adapted 2D mesh |
| `nz` | 4 | Through-thickness layers |

---

## Process Parameters

| Parameter | Symbol | Value | Units |
|---|---|---|---|
| Cutting speed | v_cut | 50 | mm/min |
| Peak power density | Q0 | 2.5 × 10¹² | W/m³ |
| Beam sigma | σ | 0.5 | mm |
| Absorption depth | 1/α | 0.1 | mm |
| Time steps | Nsteps | 60 | — |

---

## Geometry

```
        z=0  ┌─────────────────────────────────┐  top surface
             │         laser travels →          │
             │    ···●···                       │
             │                                  │  Ly = 60 mm
        z=-h └─────────────────────────────────┘  bottom surface
             x=0                             x=Lx
                       Lx = 100 mm
             plate thickness h = 1 mm
```

The laser travels along `y = Ly/2` (plate centreline) from `x = margin` to `x = Lx − margin`.

---

## Output Files

All results are saved to `D:\freefem++\laser3d\` (configurable at top of script).

| File | Description |
|---|---|
| `laser3d.pvd` | ParaView time collection — open this file |
| `frame000.vtu` … `frame059.vtu` | One VTU frame per time step |
| `burnthrough.txt` | Per-step burn-through depth log |

### Fields in each VTU frame

| Field name | Type | Description |
|---|---|---|
| `Temperature_K` | P1 nodal | Absolute temperature [K] |
| `GradT_Kpm` | P0 element | Thermal gradient magnitude [K/m] |
| `HeatIndex` | P0 element | (T − T_amb)/(T_vap − T_amb), 0→1 |
| `HAZ` | P0 element | 1 where T > 300 °C |
| `MeltPool` | P0 element | 1 where T > 635 °C (liquid Al) |
| `RecastLayer` | P0 element | 1 where 635 °C < T < 2519 °C |
| `LaserFlux` | P0 element | Volumetric source Q [W/m³] |
| `MeshSize` | P0 element | Tetrahedron volume [m³] (shows refinement) |

### Burn-through log format

```
Step  xc_mm  VapDepth_mm  MeltDepth_mm  BurnThrough  Tmax_C
1     10.13  0.00         0.00          0            1842.3
...
```

`BurnThrough = 1` when the vaporisation isotherm reaches the bottom face (z = −Lz).

---

## Requirements

| Software | Version | Notes |
|---|---|---|
| FreeFEM++ | ≥ 4.10 | [freefem.org](https://freefem.org) |
| ParaView | ≥ 5.10 | For visualisation |

FreeFEM++ plugins used (loaded automatically):

```freefem
load "iovtk"   // VTK/VTU file output
load "msh3"    // 3D mesh (buildlayers, mesh3 type)
```

No external libraries required.

---

## How to Run

### Windows
```bat
FreeFem++.exe laser_cutting_3d.edp
```

### Linux / macOS
```bash
FreeFem++ laser_cutting_3d.edp
```

### Output location
Change the output directory at the top of the script:
```freefem
string outdir = "D:\\freefem++\\laser3d\\";   // Windows
// string outdir = "/home/user/laser3d/";      // Linux/macOS
```

---

## Visualisation in ParaView

1. **File → Open** → select `laser3d.pvd` → click **Apply**
2. Change colouring to `Temperature_K` → **Rescale to Data Range**
3. Press **Play** to animate the laser traversal
4. Recommended views:

| What to see | How |
|---|---|
| HAZ rings on surface | Colour by `HAZ`, view from above (+Z) |
| Melt pool depth | **Filters → Clip** along Y axis at centreline |
| Through-thickness profile | **Filters → Slice** at x = current laser position |
| Mesh refinement moving | Colour by `MeshSize`, play animation |
| Burn-through moment | Watch `burnthrough.txt` or threshold `HeatIndex > 0.99` |

---

## Known Limitations

| Limitation | Physical meaning |
|---|---|
| No material removal | True laser cutting requires a moving domain (Stefan problem). This model keeps the mesh intact — it is a **surface heating / scribing** model, not a through-cut. |
| No melt ejection | Assist gas dynamics (N₂/O₂ jet) are not modelled |
| No latent heat | Phase-change energy (fusion/vaporisation) is not included — temperatures near T_vap are overestimated |
| 2D in-plane adaptivity only | The `adaptmesh` call refines the 2D base mesh; through-thickness resolution is uniform (`nz` layers) |
| Cold start each step | `T_old` is carried via a coarse background grid; some thermal history accuracy is lost between steps |

---

## Extending the Code

### Change material
Replace the material block at the top. Ensure all temperatures (T_sol, T_liq, T_vap, T_HAZ) are updated.

### Change laser power
```freefem
real Q0 = 5.0e12;   // increase for deeper penetration
```

### Finer through-thickness resolution
```freefem
int nz = 8;   // 8 layers instead of 4 -> 0.125 mm per layer
```

### Slower/faster cutting
```freefem
real vcut = 100./60.*1e-3;   // 100 mm/min
```

### Adapt mesh more aggressively
```freefem
int  Nadapt = 4;       // more passes
real err    = 0.002;   // tighter error target
real hmin   = sig/30.; // finer minimum element
```

---

## File Structure

```
.
├── laser_cutting_3d.edp    # Main FreeFEM++ script
├── README.md               # This file
└── results/                # Created automatically on first run
    ├── laser3d.pvd         # ParaView collection file
    ├── frame000.vtu        # Time step 0
    ├── frame001.vtu        # Time step 1
    │   ...
    ├── frame059.vtu        # Time step 59
    └── burnthrough.txt     # Depth tracking log
```

---

## Citation

If you use this code in your research, please cite:

```bibtex
@software{mishra2026laser,
  author    = {Mishra, A.},
  title     = {3D Laser Surface Heating},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.20306950},
  url       = {https://doi.org/10.5281/zenodo.20306950}
}
```

**APA:**
> Mishra, A. (2026). *3D Laser Surface Heating*. Zenodo. https://doi.org/10.5281/zenodo.20306950

**Chicago:**
> Mishra, A. 2026. "3D Laser Surface Heating." Zenodo. https://doi.org/10.5281/zenodo.20306950

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20306950.svg)](https://doi.org/10.5281/zenodo.20306950)

---

## License

MIT License — free to use, modify, and distribute with attribution.
