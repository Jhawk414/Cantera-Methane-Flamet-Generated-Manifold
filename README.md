# Cantera-Methane-Flamelet-Generated-Manifold

Creates a **Flamelet-Generated Manifold (FGM)** using a progress-variable (PV) approach for a CH₄-air counter-flow diffusion flame in Cantera/Python. Code is available for both the GRI-Mech 1.2 and 3.0 mechanisms. The code accepts user inputs of mean mixture fraction Z̃, normalized progress variable C̃, and mean mixture fraction variance Z̃'' to return β-PDF filtered combustion properties (T, ρ, Yₖ for all k species).

**Status:** Mostly-complete since May 2025. No ongoing development.

This was for my graduate-level combustion class (AE 774) completed in Spring 2025.

---

## Table of Contents
- [Repository Structure](#repository-structure)
- [Dependencies](#dependencies)
- [How to Run](#how-to-run)
- [Code Walkthrough](#code-walkthrough)
  - [Phase 1 — Flamelet Simulation](#phase-1--flamelet-simulation)
  - [Phase 2 — Z-Interpolation and Intermediate Table Assembly](#phase-2--z-interpolation-and-intermediate-table-assembly)
  - [Phase 3 — PV Normalization and C-Interpolation](#phase-3--pv-normalization-and-c-interpolation)
  - [Phase 4 — 2D Lookup Table and User Query](#phase-4--2d-lookup-table-and-user-query)
  - [Phase 5 — β-PDF Filtering and Output](#phase-5--β-pdf-filtering-and-output)
- [Notes](#notes)
- [Recommendations / Future Work](#recommendations--future-work)

---

## Repository Structure

```
AE774_CanteraManifold/
├── GRI-Mech 1.2/
│   ├── grimech12.yaml           # 32-species GRI-Mech 1.2 mechanism file
│   └── CH4 manifold gri12.py   # FGM code (32 species)
└── GRI-Mech 3.0/
    ├── gri30.yaml               # 53-species GRI-Mech 3.0 mechanism file
    └── CH4 manifold gri30.py   # FGM code (53 species)
```

The two `.py` files are functionally identical. The only differences are the `.yaml` mechanism file loaded (`grimech12.yaml` vs. `gri30.yaml`) and the species list (32 vs. 53 species).

---

## Dependencies

| Package | Purpose |
|---|---|
| `cantera` | Flamelet simulation and thermochemical property extraction |
| `numpy` | Array construction, linear spacing, transposition |
| `scipy.interpolate.interp1d` | 1D piecewise interpolation onto rectilinear grids |
| `scipy.interpolate.RegularGridInterpolator` | 2D rectilinear lookup table querying |
| `scipy.stats.beta` | β-PDF evaluation |
| `scipy.integrate.simpson` | Numerical integration of the β-PDF via composite Simpson's rule |

Install via:
```bash
pip install cantera numpy scipy
```

---

## How to Run

1. Ensure the `.yaml` mechanism file is in the same directory as the `.py` script (already the case in this repo).
2. Run from a terminal:
   ```bash
   python "CH4 manifold gri30.py"
   ```
   or for GRI 1.2:
   ```bash
   python "CH4 manifold gri12.py"
   ```
3. The script simulates and assembles the manifold — this takes several minutes. Console updates are printed as each flamelet is solved.
4. Once complete, you are prompted for three space-separated inputs:
   ```
   Z  C  Z": 0.1 0.675 0.05
   ```
   Where:
   - `Z` is the mean mixture fraction &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; (0.001 ≤ Z ≤ 0.999)
   - `C` is the normalized progress variable &nbsp; (0.001 ≤ C ≤ 0.999)
   - `Z"` is the mixture fraction variance &nbsp;&nbsp;&nbsp;&nbsp; (0.001 ≤ Z" ≤ 0.1)
5. β-PDF filtered T, ρ, and all Yₖ are printed to the console.

---

## Code Walkthrough

The code operates in five sequential phases. The output of each phase feeds directly into the next.

---

### Phase 1 — Flamelet Simulation

Eight 1D steady-state counter-flow diffusion flamelets are solved using Cantera's `CounterflowDiffusionFlame`, one per pressure in:

```python
manifold_range = [0.4, 0.5, 1, 3, 5, 7, 9, 10]   # [atm]
```

Pressure was chosen as the independent variable to span a broad progress variable solution space. Each flamelet is configured with:

| Parameter | Value |
|---|---|
| Fuel | CH₄ (pure), 300 K, ṁ = 0.15 kg/(m²·s) |
| Oxidizer | O₂:N₂ = 1:3.76 (molar), 300 K, ṁ = 0.75 kg/(m²·s) |
| Domain width | 2 cm |
| Transport model | Mixture-averaged |
| Radiation | Disabled |
| Grid refinement | ratio = 2.25, slope = 0.125, curve = 0.2, prune = 0.04 |

After each solve, Cantera returns the solution on a **non-uniform** spatial grid (more densely populated near stoichiometric conditions). The following quantities are extracted at each grid point:

| Code Variable | Symbol | Description |
|---|---|---|
| `Mix_Fraction` | Z | Mixture fraction, computed using carbon as the reference element |
| `PV` | C | Progress variable: Y\_CO₂ + Y\_H₂O |
| `Temp` | T | Temperature [K] |
| `density` | ρ | Density [kg/m³] |
| `flamelet.Y[k]` | Yₖ | Mass fraction of each of the k species |

---

### Phase 2 — Z-Interpolation and Intermediate Table Assembly

Because Cantera's output grid is non-uniform in Z, the flamelet solution must first be remapped onto an evenly-spaced mixture fraction grid before the manifold can be assembled:

```python
new_Z = np.linspace(0.01, 0.999, 75)   # 75 rectilinear Z points
```

`scipy.interpolate.interp1d` establishes the continuous functional relationship Ψ(Z) for each property Ψ ∈ {T, ρ, C, Y₁…Yₖ}, then evaluates it at every point in `new_Z`. The Z-interpolated Ψ and C for each flamelet are paired into a 75×2 sub-array `[Ψ | C]` and appended column-by-column to the growing solution matrix over the course of the pressure loop.

After all 8 flamelet iterations, this produces an intermediate solution table for each property (`ZCT_matrix`, `ZCrho_matrix`, and one dictionary entry per species in `ZCY_holder`):

**Table 1 — Z-Mapped Intermediate Solution Array**
*Shape: 75 rows × 17 columns &nbsp;(1 Z column + 8 flamelets × 2 columns each)*

| Z | Ψ (P₁=0.4 atm) | C (P₁) | Ψ (P₂=0.5 atm) | C (P₂) | ··· | Ψ (P₈=10 atm) | C (P₈) |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 0.01 | | | | | | | |
| 0.02 | | | | | | | |
| ··· | | | | | | | |
| 0.98 | | | | | | | |
| 0.99 | | | | | | | |

> Each pair of **[Ψ | C]** columns holds the Z-interpolated thermochemical property and progress variable for one flamelet at its respective pressure. The Z column is shared across all flamelets.

---

### Phase 3 — PV Normalization and C-Interpolation

At this stage, the solution is on a rectilinear Z grid, but the C axis is still non-uniform — its magnitude and range vary with each flamelet (pressure) and each Z value.

**Why normalization is needed:** Near the domain edges (min(Z) and max(Z)), the raw PV can be very small (e.g., 10⁻⁵ to 10⁻³), while near stoichiometric conditions it can exceed 0.25. Extrapolating Ψ across such a dynamic range on a common absolute-C grid would produce large errors or non-physical negative values. Normalizing C within each Z slice keeps all interpolation within the solved solution space.

The code loops over each Z "slice" (each row `zee` of the Table 1 arrays) and performs the following steps:

1. **Reshape:** Extract the row's 16 flamelet values and reshape from 1D (16 elements) into a (8 flamelets × 2 columns) array:
   ```python
   Soln_2D_T = ZCT_matrix[zee, 1:].reshape(len(manifold_range), 2)
   # Column 0: T at all 8 pressures for this Z
   # Column 1: C at all 8 pressures for this Z
   ```

2. **Normalize C** within this Z slice:
   ```python
   PV_T_norm = (Soln_2D_T[:, 1] - min(Soln_2D_T[:, 1])) \
             / (max(Soln_2D_T[:, 1]) - min(Soln_2D_T[:, 1]))
   # C_norm = (C - C_min) / (C_max - C_min)
   ```

3. **Interpolate Ψ** onto a rectilinear normalized-C grid:
   ```python
   new_C = np.linspace(0.01, 0.99, 25)    # 25 rectilinear C_norm points
   T_rectilinear = c_approx_T(new_C)      # shape: (25,)
   ```

The 25-element result for each Z slice is stacked column-by-column as the loop progresses, growing the pre-transpose lookup arrays for each property (`T_Lookup`, `rho_Lookup`, and each entry in `Y_lookup_hold`).

---

### Phase 4 — 2D Lookup Table and User Query

After looping over all 75 Z slices, each property array has shape **(25 rows × 75 columns)**. The arrays are **transposed** so that rows index Z and columns index C_norm, matching the convention expected by `RegularGridInterpolator`:

```python
T_Lookup = np.transpose(T_Lookup)   # shape: (75, 25)
```

The resulting 2D lookup table has the form:

**Table 2 — Final 2D Rectilinear Lookup Table**
*Shape: 75 rows (Z) × 25 columns (C_norm) &nbsp;— one table per property*

| Z ↓ &nbsp;\&nbsp; C → | 0.01 | 0.02 | 0.04 | ··· | 0.98 | 0.99 |
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **0.01** | | | | | | |
| **0.02** | | | | | | |
| **0.04** | | | | | | |
| **···** | | | | | | |
| **0.98** | | | | | | |
| **0.99** | | | | | | |

`RegularGridInterpolator` is built using `(new_Z, new_C)` as axis coordinates and the lookup array as data. When the user provides inputs, the interpolator performs bilinear interpolation over this 2D grid to return the unfiltered scalar property Ψ at the query point:

```python
T_reg_interp = RegularGridInterpolator((new_Z, new_C), T_Lookup)
T_unfiltered = T_reg_interp((Z_usr, C_usr))   # scalar
```

---

### Phase 5 — β-PDF Filtering and Output

The unfiltered Ψ represents a deterministic laminar-flamelet value. In turbulent combustion, the local mixture fraction fluctuates around the mean Z̃, so the FGM lookup must be weighted by a probability distribution over Z. A beta Probability Density Function (β-PDF) is used:

$$\tilde{\Psi}_i\!\left(\tilde{Z},\,\tilde{Z}''^2,\,\tilde{C}\right) = \psi_i\!\left(\tilde{Z},\tilde{C}\right) \cdot \int_0^1 \tilde{P}\!\left(Z,\tilde{Z}''^2\right) dZ$$

where the β-PDF shape is:

$$\tilde{P}\!\left(Z,\tilde{Z}''^2\right) = Z^{\alpha-1}(1-Z)^{\beta-1} \cdot \frac{\Gamma(\alpha+\beta)}{\Gamma(\alpha)\,\Gamma(\beta)}; \quad \alpha = \tilde{Z}\gamma,\;\; \beta = (1-\tilde{Z})\gamma,\;\; \gamma = \frac{\tilde{Z}(1-\tilde{Z})}{\tilde{Z}''^2} - 1$$

In code, the β-PDF is evaluated over the full `new_Z` domain using `scipy.stats.beta.pdf`. NaN values (which arise near Z = 0 or Z = 1 where shape parameters become degenerate) are removed before numerical integration with `scipy.integrate.simpson`. The resulting normalized, integrated PDF sum is a scalar correction factor multiplied against each unfiltered property:

```python
gma         = ((new_Z * (1 - new_Z)) / Zvar_usr_input**2) - 1
a           = new_Z * gma                          # α shape parameter
b           = (1 - new_Z) * gma                    # β shape parameter
bPDF        = beta.pdf(new_Z, a, b)                # β-PDF over Z domain
bPDF_no_NaN = bPDF[~np.isnan(bPDF)]               # remove boundary NaNs
P_int       = integrate.simpson(bPDF_no_NaN)       # numeric integral
bPDF_norm   = sum(bPDF_no_NaN / P_int)             # scalar correction factor

T_filtered  = T_unfiltered * bPDF_norm
```

Filtered T, ρ, and all Yₖ are printed to the console as the final output.

---

## Notes

1. A full *filtered* manifold (lookup table) of properties is not built. Instead, an "unfiltered" manifold is constructed from which properties are interpolated at the user-supplied point and then β-PDF filtered on the fly.
2. The progress variable C is normalized (C_norm = (C − C_min) / (C_max − C_min)) per Z slice so that the manifold can be built on a regular grid. This is necessary because each flamelet has a different absolute PV range depending on pressure.

---

## Recommendations / Future Work

1. **Vary flame strain rate instead of pressure.** This allows simulation of near-extinguished flamelets through fully-combusted flamelets, producing a wider and more physically meaningful manifold condition space.
2. **Filter the manifold completely.** Currently, only the properties *looked up* from the manifold are β-PDF filtered based on the user-input Z''. A fully filtered manifold would pre-integrate all property tables across Z before storage, so that filter-ready values reside directly in the lookup table.
3. **Enforce user-input bounds.** Input values are not currently validated; out-of-range inputs can cause NaNs in the β-PDF equations and extrapolation errors in `RegularGridInterpolator`. Python `exceptions` or bounds-checking guards should be added.
4. **Provide a C_min / C_max vs. Z relation as output.** Because C is normalized per Z slice, a CFD solver consuming this manifold would need the absolute C range at a given Z to convert back from a normalized progress variable to an absolute value. A supplementary lookup table or function relating C_min(Z) and C_max(Z) to Z should be added.
5. **Improve variable naming and loop consolidation.** Several intermediate variables use similar naming patterns across iterations (`holder`, `approx`, `rect` suffixes), which can create confusion when reading later code sections. A review of variable names and further consolidation of repeated operations into existing `for` loops is suggested.
