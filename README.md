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
AT2 degradation: g(d) = (1 − φ)² + k

## Reference
"T.-T. Nguyen, J. Yvonnet, D. Waldmann, and Q.-C. He. Implementation of a new strain split to model unilateral contact within the phase field method". 
International Journal for Numerical Methods in Engineering, 121(21):4717–4733, 2020. doi: https://doi.org/10.1002/nme.6463. 
URL https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.6463.

"Q.-C. He. Three-dimensional strain and stress orthogonal decompositions via an elastic energy preserving transformation." 
International Journal of Solids and Structures, 295: 112818, 2024. ISSN 0020-7683. doi: https://doi.org/10.1016/j.ijsolstr.2024.112818. 
URL https://www.sciencedirect.com/science/article/pii/S002076832400177X.
