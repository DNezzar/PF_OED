# PF_OED — Phase‑field fracture (2D Q4)

Minimal 2D phase‑field brittle fracture demo with an irreversible history field and a triangular load path: 0 → +U → 0 → −U → 0. Compares OED split vs a classic split.
[results.pdf](https://github.com/user-attachments/files/22138920/results.pdf)
<img width="800" height="500" alt="results" src="https://github.com/user-attachments/assets/88938835-be3c-4b5e-813e-4fabbb9e1d44" />

## Modes
- **SD3**: orthogonal split in the energy norm using C^(1/2) and C^(-1/2)  
  Driving energy: Ψ⁺ = ½ ε⁺ᵀ C ε⁺  
  Stress: σ = g(d) C:ε⁺ + C:ε⁻
- **Classic**: positive principal‑strain energy for damage driving  
  Stress: σ = g(d) C:ε

## Irreversibility
History variable at Gauss points: Hⁿ⁺¹(x) = max(Hⁿ(x), Ψ⁺(x))  
Nodal projection enforces φⁿ⁺¹ ≥ φⁿ

## Run
```bash
python PF_SD3.py
```

## Notes
Voigt vector: [εxx εyy γxy] with γxy = 2 εxy  
Degradation: g(d) = (1 − φ)² + k

## Reference
Implementation of a new strain split to model unilateral contact within the phase field method (IJNME 2020).
