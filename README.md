# 3D Laser Surface Heating Simulation on Al 7075-T6

<p align="center">
  <img src="https://img.shields.io/badge/FreeFEM++-Simulation-blue?style=for-the-badge&logo=gnu&logoColor=white"/>
  <img src="https://img.shields.io/badge/Al%207075--T6-Material-orange?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/3D-Adaptive%20Mesh-red?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/Beer--Lambert-Laser%20Source-purple?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/ParaView-VTK%20Export-green?style=for-the-badge"/>
  <img src="https://img.shields.io/badge/License-MIT-lightgrey?style=for-the-badge"/>
</p>

<p align="center">
  A 3D finite element simulation of <b>high-power laser surface heating</b> on a 1 mm thick
  Al 7075-T6 aerospace sheet using FreeFEM++. The laser travels along the plate centreline
  and the mesh is <b>rebuilt and adaptively refined around the moving laser spot at every
  time step</b>. A Beer-Lambert volumetric source models energy deposition through the
  plate thickness. Results are exported as a <i>dynamic PVD animation</i> for ParaView
  visualisation with full thermal, HAZ, melt pool, and recast layer fields.
</p>

![HAZ field showing adaptive mesh refinement around the moving laser spot](preview.png)

---

## Physics

The simulation solves a transient 3D heat conduction problem at each laser position,
capturing the full through-thickness thermal response. Key features include:

- Transient heat conduction PDE (backward Euler, implicit) with moving volumetric laser source
- Gaussian beam profile in the x-y plane combined with Beer-Lambert exponential decay in z
- Adaptive 2D mesh driven by a Gaussian refinement indicator at the current laser position
- Mesh extruded to 3D tetrahedra via `buildlayers` at every step — mesh follows the laser
- Previous-step temperature carried on a fixed background grid to avoid cross-mesh NaN failures
- HAZ, melt pool, recast layer, and vaporisation front tracked per element each step
- Burn-through detection: vaporisation isotherm probed along the centreline through thickness
- Real SI units throughout: metres, Watts, Kelvin, seconds

---

## Geometry

```
  3D plate (side view, x = travel direction, z = through-thickness):

  z=0   ┌─────────────────────────────────────────────┐  top surface
        │       laser travels →                        │
        │   ···●···  (Gaussian beam, sigma=0.5mm)      │
        │       ↓ Beer-Lambert decay                   │  Ly = 60 mm
        │                                              │
  z=-h  └─────────────────────────────────────────────┘  bottom surface
        x=0                                        x=Lx
                       Lx = 100 mm
        plate thickness h = 1 mm

  Top view (x-y plane):

  y=Ly  +─────────────────────────────────────────────+
        │                                             │
        │  ···●···  laser spot at y=Ly/2              │
        │                                             │
  y=0   +─────────────────────────────────────────────+
        x=0       margin=10mm          x=100mm
```

- Boundary labels 1–4 = four side faces (Dirichlet: T = T_amb)
- Top face (z = 0): laser enters via volumetric source term
- Bottom face (z = −h): heat conducted through thickness, isothermal at edges

---

## Material Parameters — Al 7075-T6

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Thermal conductivity | kc | 130 | W/(m·K) |
| Density | rho | 2810 | kg/m³ |
| Specific heat | cp | 960 | J/(kg·K) |
| Ambient temperature | Tamb | 298 (25) | K (°C) |
| Solidus | Tsol | 750 (477) | K (°C) |
| Liquidus | Tliq | 908 (635) | K (°C) |
| Vaporisation | Tvap | 2792 (2519) | K (°C) |
| HAZ threshold | THAZ | 573 (300) | K (°C) |

---

## Laser Parameters

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Peak power density | Q0 | 2.5 × 10¹² | W/m³ |
| Beam sigma | sig | 0.5 | mm |
| Absorption depth | 1/alphaAbs | 0.1 | mm |
| Cutting speed | vcut | 50 | mm/min |

---

## Process Parameters

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Plate length | Lx | 100 | mm |
| Plate width | Ly | 60 | mm |
| Plate thickness | Lz | 1 | mm |
| Start/end margin | margin | 10 | mm |
| Time steps | Nsteps | 60 | — |
| Adapt passes per step | Nadapt | 3 | — |
| Through-thickness layers | nz | 4 | — |
| Total animation frames | — | 60 | — |

---

## Governing Equations

### Transient Heat Conduction PDE

```
rho * cp * dT/dt  -  div( kc * grad(T) )  =  Q(x, y, z, t)    in Omega
```

Discretised by **backward Euler**:

```
rho*cp/dt * (T^n - T^{n-1})  -  div( kc * grad(T^n) )  =  Q^n
```

### Laser Source — Gaussian × Beer-Lambert

```
Q(x,y,z,t) = Q0 * exp( -[(x-xc)^2 + (y-yc)^2] / (2*sig^2) )
                * exp( alphaAbs * z )

z in [-Lz, 0]:
  z = 0    ->  exp(0) = 1.0   (full intensity at top surface)
  z = -Lz  ->  exp(-alphaAbs*Lz)  (attenuated at bottom face)
```

### Boundary Conditions

```
T = Tamb     on labels 1, 2, 3, 4   (all four side faces, large plate)
```

### Heat Index

```
HeatIndex = (T - Tamb) / (Tvap - Tamb)

  0.0  ->  ambient (298 K, 25 C)
  1.0  ->  vaporisation front (2792 K, 2519 C)
```

---

## Adaptive Mesh Strategy

The mesh is rebuilt and refined **every time step** around the moving laser spot:

```
Each step N:
  1. Build coarse 2D background mesh (60x40 elements over plate top face)
  2. adaptmesh x Nadapt passes using Gaussian indicator centred at (xc, yc)
       -> fine elements where the beam is intense
       -> coarse elements far from the spot
  3. buildlayers -> extrude adapted 2D mesh to 3D tetrahedra (nz layers)
  4. Interpolate T_old from fixed background grid onto new mesh
       -> clamp any NaN values to Tamb (robustness guard)
  5. Solve heat equation on adapted 3D mesh (UMFPACK direct solver)
  6. Save solution back to background grid for next step
```

### Mesh parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `Nadapt` | 3 | Adaptation passes per step |
| `err` | 0.004 | Interpolation error target (tighter = finer) |
| `hmin` | sig/20 = 0.025 mm | Finest element near beam centre |
| `hmax` | Lx/25 = 4 mm | Coarsest element far from beam |
| `nbvx` | 120 000 | Max vertices in adapted 2D mesh |
| `nz` | 4 | Through-thickness layers |

---

## Numerical Method

| Aspect | Choice |
|--------|--------|
| Spatial discretisation | Finite Element Method (FEM) |
| Element type | P1 (linear nodal, 3D) |
| Time integration | Backward Euler (implicit) |
| Linear solver | UMFPACK (direct sparse) |
| Adaptive refinement | `adaptmesh`, err=0.004, nbvx=120000 |
| Mesh per step | Rebuilt from scratch + adapted each step |
| Adaptivity iterations | 3 per step (Nadapt=3) |
| Told storage | Fixed coarse background grid (NaN-safe) |
| NaN guard | `optimize=0` on `int3d` + node-level clamp |

---

## Output Fields

Each `.vtu` frame contains the following fields:

| Field | Type | Description | Units |
|-------|------|-------------|-------|
| `Temperature_K` | P1 nodal | Absolute temperature | K |
| `GradT_Kpm` | P0 element | Thermal gradient magnitude | K/m |
| `HeatIndex` | P0 element | (T−Tamb)/(Tvap−Tamb), 0=cold 1=vap front | — |
| `HAZ` | P0 element | 1 where T > 300 °C | — |
| `MeltPool` | P0 element | 1 where T > 635 °C (liquid Al) | — |
| `RecastLayer` | P0 element | 1 where 635 °C < T < 2519 °C | — |
| `LaserFlux` | P0 element | Volumetric source Q(x,y,z) | W/m³ |
| `MeshSize` | P0 element | Tetrahedron measure (shows refinement) | m³ |

---

## Repository Structure

```
laser3d/
|
|-- laser_cutting_3d.edp          # Main FreeFEM++ simulation script
|-- README.md                     # This file
|
|-- D:\freefem++\laser3d\         # Output directory (auto-created)
    |-- laser3d.pvd               # Master animation (open in ParaView)
    |-- frame000.vtu              # Frame 0  (step 1)
    |-- frame001.vtu              # Frame 1
    |-- ...
    |-- frame059.vtu              # Frame 59 (step 60)
    |-- burnthrough.txt           # Per-step burn-through depth log
```

---

## How to Run

### Requirements

- FreeFEM++ v4.10 or later: https://freefem.org
- ParaView v5.x or later: https://www.paraview.org

### Step 1 — Run the simulation

```bash
FreeFem++ laser_cutting_3d.edp
```

The script will:
1. Print material and process parameters to console
2. For each of 60 steps: rebuild and adapt the 3D mesh around the laser spot
3. Solve the transient heat equation and extract all post-processing fields
4. Save one VTU per frame and append to master PVD file
5. Print per-frame progress including Tmax, vaporisation depth, and burn-through flag

Console output example:
```
=================================================
 3D Laser - Al 7075-T6  ADAPTIVE MOVING MESH
 Plate: 100x60x1 mm
 Speed: 50 mm/min
 Steps: 60  dt=1.6s
=================================================
Step 1/60  xc=10.00mm  Tmax=2418.3C  vapD=0.00mm  nv=4821
Step 2/60  xc=11.36mm  Tmax=2561.7C  vapD=0.62mm  nv=5104
Step 3/60  xc=12.71mm  Tmax=2748.1C  vapD=0.89mm  nv=5237
...
Step 18/60 xc=30.2mm   Tmax=2804.3C  vapD=1.00mm  nv=5318  *** BURN-THROUGH ***
```

### Step 2 — Open in ParaView

1. `File > Open` → navigate to `D:\freefem++\laser3d\`
2. Change Files of type to `All Files (*.*)`
3. Select `laser3d.pvd` → OK
4. Choose PVD Reader when prompted → OK
5. Click `Apply`
6. Set colour field to `Temperature_K`
7. Click `Rescale to Data Range Over All Timesteps`

### Step 3 — Visualise the thermal fields

**Option A — Temperature with phase boundaries**
```
Color by Temperature_K
Filters > Contour > Value = 750    (solidus, 477 C)
Filters > Contour > Value = 908    (liquidus, 635 C)
Filters > Contour > Value = 2792   (vaporisation front)
Press Play
-> Watch melt pool and vap front evolve as laser traverses plate
```

**Option B — HAZ evolution**
```
Color by HAZ  (binary: 0 or 1)
-> 1 = microstructure permanently altered (T > 300 C)
-> HAZ zone grows and persists behind laser as it travels
```

**Option C — Melt pool and recast layer**
```
Color by MeltPool   -> liquid Al region (T > 635 C)
Color by RecastLayer -> re-solidified melt (635 C < T < 2519 C)
-> Recast layer appears on kerf walls where melt re-solidifies
```

**Option D — Through-thickness profile**
```
Filters > Slice > Normal = (0,1,0), Origin = (xc, 0.030, -0.0005)
Color by Temperature_K
-> Cross-section at cut centreline shows Beer-Lambert depth gradient
-> Melt pool depth visible directly
```

**Option E — Adaptive mesh quality**
```
Color by MeshSize
-> Fine (dark blue) elements concentrated at laser spot
-> Watch the fine-mesh zone travel with the laser each frame
```

Press `Play` to watch all 60 frames of laser traversal with moving adaptive mesh.

---

## What to Look for in Results

### Beer-Lambert Depth Profile
The laser intensity decays exponentially with depth (absorption depth = 0.1 mm for Al
at 1 µm wavelength). The `Temperature_K` cross-section slice shows peak temperature
at the surface dropping to near-ambient within ~0.3–0.5 mm depth. This is physically
realistic for a focused fibre laser on aluminium.

### Adaptive Mesh Following the Laser
The `MeshSize` field visualises the mesh — finest elements (dark blue, ~0.025 mm)
are always at the current laser position, coarsening to ~4 mm at the plate edges.
In ParaView animation this looks like a spotlight of fine mesh tracking the beam.

### HAZ Growth Behind the Laser
The `HAZ` field accumulates behind the laser as it travels. Earlier positions that
were heated above 300 °C remain flagged. The HAZ width reflects the beam sigma
and travel speed — faster cutting = narrower HAZ.

### Burn-Through Moment
When the vaporisation isotherm reaches z = −1 mm (bottom face), the console prints
`*** BURN-THROUGH ***` and `burnthrough.txt` records `BurnThrough=1`. In ParaView,
threshold on `HeatIndex > 0.99` to see the vap-front isosurface approaching the
bottom face.

---

## Common Errors and Fixes

| Error | Cause | Fix |
|-------|-------|-----|
| `syntax error before token _` | Underscore in identifier | Use camelCase: `alphaAbs` not `alpha_abs` |
| `nan != nan diff nan` | Cross-mesh interpolation returns NaN | `optimize=0` on `int3d` + node-level NaN clamp |
| `Assertion fail nrfmid` | `buildlayers` label arrays wrong size | Remove `labelmid/labelup/labeldown`; use simple 4-arg form |
| `No operator .volume` | `.volume` not valid on `mesh3` elements | Use `.measure` instead |
| All frames same temperature | `Told` not updating | Ensure `Told = ui` executes inside the time loop |
| Mesh too coarse near spot | `err` too large or `nbvx` too low | Set `err=0.004`, `nbvx=120000`, `hmin=sig/20` |
| Non-ASCII compile error | UTF-8 dash or accented char in string | Use only 7-bit ASCII in all string literals |

---

## Extending the Model

| Extension | What to change |
|-----------|----------------|
| Different material | Update kc, rho, cp, Tsol, Tliq, Tvap, THAZ |
| Higher power / deeper cut | Increase Q0; raise alphaAbs for shallower absorption |
| Faster travel | Increase vcut; dt updates automatically |
| Finer mesh near spot | Decrease err to 0.002, hmin to sig/30 |
| More through-thickness resolution | Increase nz to 8 (0.125 mm/layer) |
| Persistent thermal history | Use a finer background grid (increase nxB, nyB, nzB) |
| Convective surface cooling | Add Robin BC on top/bottom faces via `int2d` |
| Curved cut path | Parameterise xc, yc as functions of step index |

---

## Citation

If you use this code in your research, please cite:

```bibtex
@software{mishra2026laser3d,
  author    = {Mishra, A.},
  title     = {3D Laser Surface Heating},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.20306950},
  url       = {https://doi.org/10.5281/zenodo.20306950}
}
```

Plain text citation:

> Mishra, A. (2026). *3D Laser Surface Heating*. Zenodo. https://doi.org/10.5281/zenodo.20306950

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20306950.svg)](https://doi.org/10.5281/zenodo.20306950)

---

## Author

**akshansh11**
GitHub: https://github.com/akshansh11

---

## License

<p>
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">
<img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png"/>
</a>
<br/>
This work is licensed under a
<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.
</p>

You are free to:

- **Share** — copy and redistribute the material in any medium or format
- **Adapt** — remix, transform, and build upon the material

Under the following terms:

- **Attribution** — You must give appropriate credit to akshansh11 and provide a link to this repository
- **NonCommercial** — You may not use the material for commercial purposes

Copyright 2026 akshansh11. All rights reserved for commercial use.
