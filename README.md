# PF_OED — Phase-field brittle fracture with Orthogonal Energy Decomposition

Minimal 2D (Q4) phase-field demo comparing **OED** (a.k.a. SD3) vs a **classic** tensile split.  
It includes a monotone **history variable** H for irreversibility, a **nodal projection** (dⁿ⁺¹ ≥ dⁿ), and a **triangular loading**: 0 → +U → 0 → −U → 0.

---

## What’s inside

- **OED split (energy-orthogonal)** with two options:
  - `oed_method="projector"`: closed-form 2×2 (invariants + projectors), no eigvecs
  - `oed_method="eig"`: standard eigendecomposition in the transformed space
- **Classic split** (reference): positive principal strains of ε (no C¹ᐟ²) and degradation on total stress
- **Plane strain**, **isotropic elasticity** in Voigt `[εxx, εyy, γxy]` with `γxy = 2 εxy`
- **Unilateral stress** for OED:  
  `σ = g(d) C ε⁺ + C ε⁻`, with `g(d) = (1 − d)² + k`
- **Outputs**: figures + `.dat` files for post-processing

> This is a didactic single-element example (Q4, 4 Gauss points). All model parameters are grouped at the **top of the script**.

---

## Quick start

```bash
# minimal environment
python -m venv .venv
source .venv/bin/activate  # (Windows: .venv\Scripts\activate)
pip install numpy matplotlib
python PF_OED.py
