**⭐ If you find this code useful for your research or interesting, please consider giving it a star! and share with your friends**

**🎓 This code is not just for research - it's also fascinating to watch dendrites grow in real-time!**
**Whether you're a materials science student, a physics enthusiast, or just curious about pattern formation, watching six-fold symmetric dendrites grow from a tiny seed is mesmerizing.**

# Concentration-Phase Field Method for Dendritic Growth

An educational Concentration-phase-field method implementation perfect for PFM beginners to learn detailed simulation techniques. Based on the anisotropic Allen-Cahn equation, solved using finite difference method (forward Euler in time, five-point central difference in space), implementing dimensionless phase-field method with rigorous concentration-field coupling for six-fold / four-fold symmetric dendritic growth in binary alloy solidification.

**Key Features for Learning**:
- Clear finite difference discretization
- Explicit time integration scheme
- Modular code structure for easy understanding

## Example Videos and Images

### Example 1: Basic Dendritic Growth

![Basic Growth Example](examples/images/basic_growth_preview.png)
*Video: Combined view of phase, concentration, and C/C∞ fields*

<video src="https://github.com/user-attachments/assets/6628695f-00f5-48e8-8ad1-93370805df96" controls="controls" width="100%"></video>

Six-fold symmetric dendrite growing from a central seed. The video shows simultaneous evolution of three fields in subplots.

**Performance**: GPU parallel computing completes the simulation in just **391ms**, achieving a **388× speedup** compared to MATLAB CPU single-core serial computing (152s).

**Test Environment**:
- **CPU**: AMD Ryzen 7 8845HS w/ Radeon 780M Graphics @ 3.80 GHz
- **GPU**: 4GB NVIDIA GeForce RTX 3050 Laptop GPU

### Example 2: Velocity Field Effects

![Velocity Field Example](examples/images/velocity_field_preview.png)
*Videos: Combined fields, Streamlines with velocity magnitude background*

<video src="https://github.com/user-attachments/assets/6e780264-69dc-40a3-af7b-bf80946dcc0e" controls="controls" width="100%"></video>

<video src="https://github.com/user-attachments/assets/77b1f14a-9711-406f-befa-d10c1243ebac" controls="controls" width="100%"></video>

Dendrite growing under uniform flow (forced convection). The flow breaks the six-fold symmetry and creates asymmetric growth patterns:
- **Upstream side**: Fresh fluid brings lower solute concentration → **Faster growth**
- **Downstream side**: Solute accumulates → Higher concentration → **Slower growth**

### Parameter Dependence

![Parameter Dependence](examples/images/parameter_dependence.png)
*Tip curvature radius and growth velocity as functions of material parameters*

Tip curvature radius and growth velocity vary with λ

---

## 🚀 CUDA GPU Demonstration

**Million-grid computing**: 1000×1000 grid, 50τ₀ simulation completed in just **22.3 seconds**

<video src="https://github.com/user-attachments/assets/d53e2fdb-9c83-46bc-b3e9-4b438198939c" controls="controls" width="100%"></video>

**Key Features**:
- **Massive speedup**: 388× faster than MATLAB (391ms vs 152s for basic simulation)
- **Complete physics**: Phase field (5 terms) + Concentration field (3 terms)
- **High performance**: Optimized for RTX 30-series GPUs
- **Easy visualization**: MATLAB integration for results analysis

**Quick Start**:
```powershell
cd CUDA
.\build_ninja_devcmd.bat
```

⚠️ **Note**: CUDA version requires complex environment configuration (CUDA Toolkit, CMake, compilers). If you encounter issues running the CUDA version, it is most likely due to incorrect environment configuration. **For beginners, we recommend using the MATLAB version** for easier setup and usage.

---

### Troubleshooting Guide

If you encounter errors, bugs, or unphysical phenomena:

**Step 1**: Check [Known Limitations](#known-limitations-) ⭐
- Many issues are documented here with solutions

**Step 2**: Check [Common Issues](#common-issues)
- Quick fixes for frequent problems

**Step 3**: Check [Velocity Field Physical Limitations](#velocity-field-physical-limitations-) ⭐
- Important notes about velocity field implementation

**Step 4**: Report on [GitHub Issues](https://github.com/pavedplaza/phase-field-dendrite-basic/issues)

---

## Table of Contents

- [Basics](#basics)
  - [Quick Start ⭐](#quick-start-)
  - [Common Issues](#common-issues)
  - [Project Structure](#project-structure)
  - [Theory Background](#theory-background)
  - [Default Parameters](#default-parameters)
  - [Physical Units](#physical-units)
  - [Visualization & Output Recommendations](#visualization--output-recommendations)
- [Advanced](#advanced)
  - [Video Output Configuration](#video-output-configuration)
  - [Velocity Field Physical Limitations ⭐](#velocity-field-physical-limitations-)
  - [Known Limitations ⭐](#known-limitations-)
  - [Initial Seed Configuration](#initial-seed-configuration)

---

## Basics

### Quick Start ⭐

#### 3 Steps to Run

**Step 1**: Download project
```bash
git clone https://github.com/pavedplaza/phase-field-dendrite-basic.git
cd phase-field-dendrite-basic
```

**Step 2**: Open MATLAB and add to path
```matlab
cd 'path/to/phase-field-dendrite-basic'
addpath(genpath('.'))
```

**Step 3**: Run simulation
```matlab
% Method 1: Use main script (recommended for beginners)
run('phase_field_simulation.m')

% Method 2: Use function interface (flexible configuration)
params = config_default();
run_simulation(params);
```

**✅ That's it! Wait for simulation to complete (very soon)**

Results will be saved in `simulation_data.mat` with analysis plots generated.

### Common Issues

1. **Running too slowly**: Set `params.visualization = false`
2. **Has 'enable_resampling = true' not taken effect?**: `enable_resampling` only takes effect when `tip_tracking = true`

---

### Project Structure

This section explains the main files and core functions in the project.

#### Main Programs

**`run_simulation.m`** - Main simulation function
- Core simulation loop with time integration
- Handles all physics calculations (phase field, concentration field, velocity field)
- Manages data saving and visualization
- **Usage**: `results = run_simulation(params);`

**`config_default.m`** - Default parameter configuration
- Sets all simulation parameters (material properties, numerical settings, boundary conditions)
- Customizable before running simulation
- **Usage**: `params = config_default();`

**`phase_field_simulation.m`** - Main script for beginners
- Simple entry point for running simulations
- Automatically loads default config and runs simulation
- **Usage**: `run('phase_field_simulation.m')`

#### PFM_core Functions

All functions in `PFM_core/` folder handle tip tracking and analysis:

**`find_six_tips_direction_specific.m`** - Detect six dendrite tips
- Searches for tips along six directions (0°, 60°, 120°, 180°, 240°, 300°)
- Requires: `velocity_field = 'none'`, `m_aniso = 6`, central seed

**`calculate_six_tip_dynamics.m`** - Calculate tip dynamics
- Computes tip position, velocity, curvature radius
- Maintains history buffer for smoothing

**`calculate_instantaneous_tip_velocity.m`** - Instantaneous tip velocity
- Calculates velocity from tip position changes

**`resample_contour_equidistant.m`** - Resample interface contour
- Resamples contour points at equal intervals
- Required for accurate tip velocity calculation

**`calculate_contour_line_intersections.m`** - Line-contour intersections
- Calculates intersection points between search rays and contour

**`generate_tip_detailed_analysis.m`** - Generate detailed analysis
- Creates comprehensive tip tracking reports

**`plot_six_tip_dynamics_analysis.m`** - Plot tip dynamics
- Generates analysis plots for tip position, velocity, curvature

---

### Theory Background

#### Phase Field Equation

```
A(ψ)² * [ k * (1 + (1 - k) * U) ] ∂φ/∂t  = ∇·[ A(ψ)² ∇φ ]
            - ∂/∂x [ A(ψ) A'(ψ) ∂φ/∂y ]
            + ∂/∂y [ A(ψ) A'(ψ) ∂φ/∂x ]
            + φ (1 - φ²)
            - λ (1 - φ²)² ( θ + k * U )
```

Where:
- `ϕ` (phi): Phase field (ϕ=1 solid, ϕ=-1 liquid)
- `k`: Partition coefficient
- `λ` (lambda): Coupling constant
- `A(ψ)`: Anisotropy function

#### Concentration Field Equation

```
[ (1+k) - (1-k)φ ] / 2 * ∂U/∂t =
    ∇·[ D * (1-φ)/2 * ∇U
        + (1 + (1-k)U)/(2√2) * (∂φ/∂t) * (∇φ / |∇φ|) ]
    + [1 + (1-k)U]/2 * ∂φ/∂t
    - (1-φ)(1-k)/4 * v * {
        [1 + k - (1-k)φ] ∇U - [1 + (1-k)U] ∇φ
      }
```

Where:
- `U`: Dimensionless concentration
- `D`: Diffusion coefficient
- `v`: Velocity field

---

#### Numerical Method

- **Spatial discretization**: Finite difference with 5-point stencil
- **Time integration**: Explicit Euler scheme  
- **Time step**: Calculated based on CFL condition: `dt = safety_factor * dx² / (4D)`
- **Boundary conditions**: Periodic or zero-gradient

### Default Parameters

#### Material Parameters

| Parameter | Symbol | Default | Description |
|-----------|--------|---------|-------------|
| Partition coefficient | k | 0.15 | Solute partition coefficient |
| Anisotropy strength | ε | 0.02 | Strength of anisotropy |
| Anisotropy mode | m | 6 | Mode of anisotropy (6-fold) |
| Coupling constant | λ | 10.0 | Coupling strength |
| Diffusion coefficient | D | 6.267 | Dimensionless diffusivity |
| Undercooling | θ | -0.2 | Dimensionless undercooling |

#### Numerical Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| Grid size (Nx, Ny) | 300×300 | Grid resolution |
| Spatial step (dx, dy) | 0.8 W | Grid spacing |
| Total time | 10.0 τ₀ | Simulation duration |
| Output interval | 0.5 τ₀ | Data save frequency |

---

### Physical Units

**Important**: This simulation uses **dimensionless phase field method**. All units are normalized:

- **Time unit**: `τ₀` (relaxation time)
- **Length unit**: `W` (interface width)
- `dt` (time step): Default ≈ 0.0064 τ₀ (calculated from CFL condition)
- `dx` (grid spacing): Default = 0.8 W
- Seed crystal radius: Default = 3.2 W

---

### Visualization & Output Recommendations

#### Scenario 1: Real-time Observation (Visualization Mode)

```matlab
params.visualization = true;      % Enable real-time display
params.tip_tracking = false;       % Disable tip tracking
params.output_interval = 0.5;      % Larger interval for smooth visualization

% Recommended:
% - Grid: 300×300 (faster rendering)
% - Total time: 5-10 τ₀
% - Use: Qualitative observation, demonstrations
```

#### Scenario 2: Data Analysis (Non-visualization Mode)

```matlab
params.visualization = false;     % Disable visualization for speed
params.tip_tracking = true;       % Enable tip tracking
params.output_interval = 0.01;    % Small interval for high-resolution data

% Recommended:
% - Grid: 300×400 or larger
% - Total time: 10-20 τ₀
% - Use: Quantitative analysis, tip velocity measurement, publication
```

**Why?**
- **Tip tracking requires high temporal resolution**: Tip velocity changes rapidly
- **Visualization is computationally expensive**: Rendering slows down computation
- **Post-processing plots are always generated**: Even with `visualization=false`

---

## Advanced

### Video Output Configuration

The simulation can generate MP4 videos showing the evolution of multiple fields.

#### How to Enable Video Generation

```matlab
params = config_default();
params.save_video = true;          % Enable video saving
params.video_fps = 10;             % Video frames per second (default: 10)
params.visualization = true;       % Required: must enable visualization
params.output_interval = 0.5;      % Controls video smoothness

results = run_simulation(params);
```

**Important**:
- `save_video` only works when `visualization = true`
- Videos are generated after simulation completes (post-processing)
- Video generation does NOT slow down the simulation itself

#### Output Videos

**For basic growth** (no velocity field):
- `basic_growth_combined.mp4`: Single video with 3 subplots showing:
  - Phase field (ϕ)
  - Concentration field (U)
  - C/C∞ field

**For velocity field simulations**:
- `velocity_field_combined.mp4`: Combined fields (same as above)
- `velocity_field_contour.mp4`: Velocity field quiver plot with dendrite contour

#### Performance Tips

- **Larger `output_interval`** → Fewer frames → Faster video generation, smoother playback
- **Smaller `output_interval`** → More frames → Higher temporal resolution, slower generation
- **Recommended**: `output_interval = 0.5` for good balance

---

### Velocity Field Physical Limitations ⭐

**Important**: The current implementation does NOT solve the full Navier-Stokes equations.

This means:

- **All flow types** (uniform, shear, vortex): Use a **simple permeability model** to prevent fluid penetration into solid phase
  - Velocity in solid phase is set to zero via `(1-φ)²` permeability factor

- **Shear flow** (`velocity_field = 'shear'`): Shows a **linear velocity gradient** without proper rotational physics

- **Vortex flow** (`velocity_field = 'vortex'`): Shows a **simple rotational velocity field** without considering:
  - Conservation of angular momentum
  - Viscous dissipation effects
  - Navier-Stokes equations

For physically accurate fluid-solid interaction simulations, a full CFD (Computational Fluid Dynamics) implementation with two-way coupling would be required.

**Recommendation**:
- Use **uniform flow** (`velocity_field = 'uniform'`) for most studies (simplest approximation)
- The shear and vortex options are provided for **educational demonstrations only**

---

### Known Limitations ⭐

#### Tip Tracking Constraints

Six tips are detected using:

1. **Contour extraction**: MATLAB `contour` function at ϕ=0
2. **Direction-specific search**: Six rays at 60° intervals
3. **Intersection calculation**: Ray-contour intersection points
4. **Validation**: Distance and angle deviation checks

**Important**: Due to algorithmic limitations, tip tracking (`params.tip_tracking`) currently has **three strict requirements**:

1. **No velocity field**: Only `velocity_field = 'none'` is supported
2. **Six-fold symmetry**: Anisotropy mode must be `m_aniso = 6`
3. **Central seed**: Seed crystal must be at domain center

```matlab
% ✅ CORRECT: Tip tracking with all requirements met
params.tip_tracking = true;
params.velocity_field = 'none';       % Only 'none' is supported
params.m_aniso = 6;                   % Must be 6 for six-fold symmetry

% ❌ INCORRECT 1: Tip tracking with velocity field (NOT SUPPORTED)
params.tip_tracking = true;
params.velocity_field = 'uniform';    % Will cause errors

% ❌ INCORRECT 2: Tip tracking with wrong anisotropy mode (NOT SUPPORTED)
params.tip_tracking = true;
params.m_aniso = 4;                   % Must be 6, not 4
```

**Why?**
- The tip detection algorithm assumes **six-fold symmetric growth from a central seed**
- **Velocity fields break this symmetry**, making tip detection unreliable
- **Seed must be at domain center** for the angle-based tip detection algorithm

---

### Initial Seed Configuration

**Fixed assumptions**:
- Seed crystal is **always initialized at the domain center**
- Growth pattern is **six-fold symmetric** (by default, `params.m_aniso = 6`)
- Seed radius: Default = 3.2 W

**To modify these assumptions**, you would need to:
1. Change the seed initialization in `run_simulation.m` (lines ~160-180)
2. Update the tip detection algorithm in `PFM_core/detect_six_tips.m`
3. Adjust target angles in `params.six_tip_angles` (default: 0:60:300 degrees)

**Recommendation**:
- For tip tracking studies: Always use `params.velocity_field = 'none'`
- For velocity field studies: Disable tip tracking (`params.tip_tracking = false`)

---

---

**Last Updated**: 2026-02-27

**Maintainer**: pavedplaza <2300837983@qq.com>

**Affiliation**: China Zhejiang University(ZJU) And Chongqing University(CQU)
