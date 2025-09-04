# -*- coding: utf-8 -*-
"""
Q4 2D - Phase-field + OED vs Classique
OED = Orthogonal Energy Decomposition (anciennement SD3)
Options OED:
  - oed_method="projector": split OED en forme fermée 2D (invariants + projecteurs)
  - oed_method="eig"      : split OED par décomposition en valeurs propres
Irréversibilité: variable historique H et projection nodale phi^{n+1} >= phi^n
Chargement triangulaire: 0 -> +U -> 0 -> -U -> 0
"""

import numpy as np
import math as mp
import matplotlib.pyplot as plt

# ===================== PARAMÈTRES (tout au début) ===================== #
# --- Maillage (Q4, 1 élément) ---
lx, ly = 1.0, 1.0
nx, ny = 1, 1
dim_n, dim_e = 2, 4

# --- Matériau (déformation plane) ---
E, nu, h = 1.0, 0.3, 1.0

# --- Phase-field ---
Gc, l0 = 1.0, 1.0

# --- Chargement triangulaire ---
# 0->+U, +U->0, 0->-U, -U->0
U_imposed = 2.0
n1, n2, n3, n4 = 25, 25, 25, 25  #Nb incréments par phase

# --- Numérique ---
k_res = 1e-8  # résiduel pour éviter rigidité nulle
# ===================================================================== #

# ------------- Maillage Q4 : 1 élément -------------
nelem = nx * ny
nnode = (nx + 1) * (ny + 1)
ncoor = np.zeros((nnode, dim_n))
connec = np.zeros((nelem, dim_e), dtype=int)

n = 0
for j in range(ny + 1):
    for i in range(nx + 1):
        ncoor[n, 0] = i * lx / nx
        ncoor[n, 1] = j * ly / ny
        n += 1
e = 0
for j in range(1, ny + 1):
    for i in range(nx):
        connec[e, 0] = i + (j - 1) * (nx + 1)
        connec[e, 1] = i + 1 + (j - 1) * (nx + 1)
        connec[e, 2] = i + 1 + j * (nx + 1)
        connec[e, 3] = i + j * (nx + 1)
        e += 1

# ------------- Matériau (déformation plane) -------------
lam = (E * nu) / ((1 + nu) * (1 - 2 * nu)) if (1 - 2*nu) != 0 else 0.0
mu_shear = E / (2 * (1 + nu))  # μ

# Voigt [xx, yy, γxy] avec γ=2*εxy
C = np.zeros((3, 3))
C[0, 0] = lam + 2 * mu_shear
C[1, 1] = lam + 2 * mu_shear
C[0, 1] = lam
C[1, 0] = lam
C[2, 2] = mu_shear

def compute_C_sqrt_matrices():
    """Matrices C^(1/2) et C^(-1/2) (formules isotropes 2D)."""
    kappa = lam + 2 * mu_shear / 3.0
    C_half = np.zeros((3, 3))
    sqrt_k = np.sqrt(max(kappa, 0.0))
    sqrt_mu = np.sqrt(max(mu_shear, 0.0))
    C_half[0, 0] = (sqrt_k + sqrt_mu) / np.sqrt(2)
    C_half[1, 1] = (sqrt_k + sqrt_mu) / np.sqrt(2)
    C_half[0, 1] = (sqrt_k - sqrt_mu) / np.sqrt(2)
    C_half[1, 0] = (sqrt_k - sqrt_mu) / np.sqrt(2)
    C_half[2, 2] = np.sqrt(2 * mu_shear)

    C_mhalf = np.zeros((3, 3))
    if sqrt_k > 0 and sqrt_mu > 0:
        C_mhalf[0, 0] = 1.0 / (2 * np.sqrt(2 * kappa)) + 1.0 / (2 * np.sqrt(2 * mu_shear))
        C_mhalf[1, 1] = 1.0 / (2 * np.sqrt(2 * kappa)) + 1.0 / (2 * np.sqrt(2 * mu_shear))
        C_mhalf[0, 1] = 1.0 / (2 * np.sqrt(2 * kappa)) - 1.0 / (2 * np.sqrt(2 * mu_shear))
        C_mhalf[1, 0] = 1.0 / (2 * np.sqrt(2 * kappa)) - 1.0 / (2 * np.sqrt(2 * mu_shear))
        C_mhalf[2, 2] = 1.0 / np.sqrt(2 * mu_shear)
    return C_half, C_mhalf

C_half, C_mhalf = compute_C_sqrt_matrices()

# ------------- B-matrix & rigidité locale -------------
def Strain(xsi, eta, x, y):
    J = np.zeros((2, 2))
    N = np.zeros((4, 8))
    I = np.zeros((3, 4))
    Jaco = np.zeros((4, 4))

    J[0, 0] = 0.25 * ((x[1] - x[0]) * (1 - eta) + (x[2] - x[3]) * (1 + eta))
    J[0, 1] = 0.25 * ((y[1] - y[0]) * (1 - eta) + (y[2] - y[3]) * (1 + eta))
    J[1, 0] = 0.25 * ((x[3] - x[0]) * (1 - xsi) + (x[2] - x[1]) * (1 + xsi))
    J[1, 1] = 0.25 * ((y[3] - y[0]) * (1 - xsi) + (y[2] - y[1]) * (1 + xsi))

    N[0, 0] = -0.25 * (1 - eta); N[0, 2] = 0.25 * (1 - eta)
    N[0, 4] =  0.25 * (1 + eta); N[0, 6] = -0.25 * (1 + eta)
    N[1, 0] = -0.25 * (1 - xsi); N[1, 2] = -0.25 * (1 + xsi)
    N[1, 4] =  0.25 * (1 + xsi); N[1, 6] =  0.25 * (1 - xsi)
    N[2, 1] = -0.25 * (1 - eta); N[2, 3] =  0.25 * (1 - eta)
    N[2, 5] =  0.25 * (1 + eta); N[2, 7] = -0.25 * (1 + eta)
    N[3, 1] = -0.25 * (1 - xsi); N[3, 3] = -0.25 * (1 + xsi)
    N[3, 5] =  0.25 * (1 + xsi); N[3, 7] =  0.25 * (1 - xsi)

    Jinv = np.linalg.inv(J)
    Jaco[0, 0] = Jinv[0, 0]; Jaco[0, 1] = Jinv[0, 1]
    Jaco[1, 0] = Jinv[1, 0]; Jaco[1, 1] = Jinv[1, 1]
    Jaco[2, 2] = Jinv[0, 0]; Jaco[2, 3] = Jinv[0, 1]
    Jaco[3, 2] = Jinv[1, 0]; Jaco[3, 3] = Jinv[1, 1]

    I[0, 0] = 1; I[1, 3] = 1; I[2, 1] = 1; I[2, 2] = 1
    B = I @ (Jaco @ N)
    return B, C

# ------------- OED: deux méthodes -------------
def oed_split_eps_eig(eps_voigt, C_half, C_mhalf):
    """OED via eigendecomposition classique."""
    eps_tilde = C_half @ eps_voigt
    E = np.array([[eps_tilde[0], eps_tilde[2] / 2.0],
                  [eps_tilde[2] / 2.0, eps_tilde[1]]])
    w, V = np.linalg.eig(E)
    idx = w.argsort()[::-1]
    w, V = w[idx], V[:, idx]
    w_plus = np.maximum(w, 0.0); w_minus = np.minimum(w, 0.0)
    Ep = np.zeros((2, 2)); Em = np.zeros((2, 2))
    for i in range(2):
        n = V[:, i].reshape(-1, 1)
        Ep += w_plus[i] * (n @ n.T)
        Em += w_minus[i] * (n @ n.T)
    v = lambda M: np.array([M[0, 0], M[1, 1], 2.0 * M[0, 1]])
    eps_plus  = C_mhalf @ v(Ep)
    eps_minus = C_mhalf @ v(Em)
    return eps_plus, eps_minus

def oed_split_eps_projector(eps_voigt, C_half, C_mhalf, tol=1e-12):
    """OED en forme fermée 2D (invariants + projecteurs, sans eigvec)."""
    eps_t = C_half @ eps_voigt
    T = np.array([[eps_t[0], eps_t[2] / 2.0],
                  [eps_t[2] / 2.0, eps_t[1]]])
    tr  = T[0, 0] + T[1, 1]
    det = T[0, 0] * T[1, 1] - T[0, 1] * T[1, 0]
    disc = tr * tr - 4.0 * det
    if disc < 0.0:
        disc = 0.0
    root = np.sqrt(disc)
    lam1 = 0.5 * (tr + root)
    lam2 = 0.5 * (tr - root)

    I2 = np.eye(2)
    if abs(lam1 - lam2) > tol:
        E1 = (T - lam2 * I2) / (lam1 - lam2)
        E2 = I2 - E1
    else:
        if lam1 >= 0.0:
            E1, E2 = I2, np.zeros((2, 2))
        else:
            E1, E2 = np.zeros((2, 2)), I2
        lam2 = lam1

    Tplus  = max(lam1, 0.0) * E1 + max(lam2, 0.0) * E2
    Tminus = min(lam1, 0.0) * E1 + min(lam2, 0.0) * E2

    v = lambda M: np.array([M[0, 0], M[1, 1], 2.0 * M[0, 1]])
    eps_plus  = C_mhalf @ v(Tplus)
    eps_minus = C_mhalf @ v(Tminus)
    return eps_plus, eps_minus

def oed_split_eps(eps_voigt, method="projector"):
    if method == "projector":
        return oed_split_eps_projector(eps_voigt, C_half, C_mhalf)
    elif method == "eig":
        return oed_split_eps_eig(eps_voigt, C_half, C_mhalf)
    else:
        raise ValueError("oed_method must be 'projector' or 'eig'.")

# ------------- Énergie positive "classique" -------------
def classic_energy_positive(eps_voigt):
    exx, eyy, gxy = eps_voigt
    exy = gxy / 2.0
    E2 = np.array([[exx, exy], [exy, eyy]])
    w, V = np.linalg.eig(E2)
    idx = w.argsort()[::-1]
    w = w[idx]; V = V[:, idx]
    w_plus = np.maximum(w, 0.0)
    Ep = np.zeros((2, 2))
    for i in range(2):
        n = V[:, i].reshape(-1, 1)
        Ep += w_plus[i] * (n @ n.T)
    eps_plus_voigt = np.array([Ep[0, 0], Ep[1, 1], 2.0 * Ep[0, 1]])
    Psi_plus = 0.5 * (eps_plus_voigt @ (C @ eps_plus_voigt))
    return Psi_plus

# ------------- Phase-field (H comme driver) -------------
def K_phi_elem(x, y, eta, xsi, H_scalar):
    N_phi = np.zeros((1, 4))
    dN_phi = np.zeros((2, 4))
    J = np.zeros((2, 2))

    N_phi[0, 0] = 0.25 * (1 - xsi) * (1 - eta)
    N_phi[0, 1] = 0.25 * (1 + xsi) * (1 - eta)
    N_phi[0, 2] = 0.25 * (1 + xsi) * (1 + eta)
    N_phi[0, 3] = 0.25 * (1 - xsi) * (1 + eta)

    J[0, 0] = 0.25 * ((x[1] - x[0]) * (1 - eta) + (x[2] - x[3]) * (1 + eta))
    J[0, 1] = 0.25 * ((y[1] - y[0]) * (1 - eta) + (y[2] - y[3]) * (1 + eta))
    J[1, 0] = 0.25 * ((x[3] - x[0]) * (1 - xsi) + (x[2] - x[1]) * (1 + xsi))
    J[1, 1] = 0.25 * ((y[3] - y[0]) * (1 - xsi) + (y[2] - y[1]) * (1 + xsi))

    dN_phi[0, 0] = -0.25 * (1 - eta); dN_phi[0, 1] =  0.25 * (1 - eta)
    dN_phi[0, 2] =  0.25 * (1 + eta); dN_phi[0, 3] = -0.25 * (1 + eta)
    dN_phi[1, 0] = -0.25 * (1 - xsi); dN_phi[1, 1] = -0.25 * (1 + xsi)
    dN_phi[1, 2] =  0.25 * (1 + xsi); dN_phi[1, 3] =  0.25 * (1 - xsi)

    detJ = np.linalg.det(J)
    B_phi = np.linalg.inv(J) @ dN_phi
    return (Gc * l0 * (B_phi.T @ B_phi) + ((Gc / l0) + 2 * H_scalar) * (N_phi.T @ N_phi)) * detJ

def K_phi_integ(x, y, H_vec):
    gp = 1 / mp.sqrt(3)
    xsi = [-gp, gp, gp, -gp]
    eta = [-gp, -gp, gp, gp]
    K = np.zeros((4, 4))
    for i in range(4):
        K += K_phi_elem(x, y, eta[i], xsi[i], H_vec[i, 0])
    return K

def Resi_elem(x, y, eta, xsi, H_scalar):
    N_phi = np.zeros((1, 4))
    J = np.zeros((2, 2))

    N_phi[0, 0] = 0.25 * (1 - xsi) * (1 - eta)
    N_phi[0, 1] = 0.25 * (1 + xsi) * (1 - eta)
    N_phi[0, 2] = 0.25 * (1 + xsi) * (1 + eta)
    N_phi[0, 3] = 0.25 * (1 - xsi) * (1 + eta)

    J[0, 0] = 0.25 * ((x[1] - x[0]) * (1 - eta) + (x[2] - x[3]) * (1 + eta))
    J[0, 1] = 0.25 * ((y[1] - y[0]) * (1 - eta) + (y[2] - y[3]) * (1 + eta))
    J[1, 0] = 0.25 * ((x[3] - x[0]) * (1 - xsi) + (x[2] - x[1]) * (1 + xsi))
    J[1, 1] = 0.25 * ((y[3] - y[0]) * (1 - xsi) + (y[2] - y[1]) * (1 + xsi))
    detJ = np.linalg.det(J)
    return (2 * H_scalar) * N_phi * detJ

def Resi_integ(x, y, H_vec):
    gp = 1 / mp.sqrt(3)
    xsi = [-gp, gp, gp, -gp]
    eta = [-gp, -gp, gp, gp]
    R = np.zeros((1, 4))
    for i in range(4):
        R += Resi_elem(x, y, eta[i], xsi[i], H_vec[i, 0])
    return R

# ------------- Chargement triangulaire -------------
u2_history = np.concatenate([
    np.linspace(0.0,  +U_imposed, n1, endpoint=False),
    np.linspace(+U_imposed, 0.0,  n2, endpoint=False),
    np.linspace(0.0,  -U_imposed, n3, endpoint=False),
    np.linspace(-U_imposed, 0.0,  n4, endpoint=True)
])
nb_inc_total = u2_history.size

# ------------- Helper pour axes à 0 -------------
def add_zero_axes(ax):
    ax.axhline(0, color='k', ls='--', lw=1, alpha=0.6)
    ax.axvline(0, color='k', ls='--', lw=1, alpha=0.6)
    xmin, xmax = ax.get_xlim(); ymin, ymax = ax.get_ylim()
    ax.set_xlim(min(xmin, 0), max(xmax, 0))
    ax.set_ylim(min(ymin, 0), max(ymax, 0))

# ------------- Simulation avec options -------------
def run_simulation(mode="oed", oed_method="projector"):
    """
    mode='oed' ou 'classic'
    oed_method='projector' (forme fermée) ou 'eig' (décomp. propre)
    """
    assert mode in ("oed", "classic")
    assert oed_method in ("projector", "eig")

    x = ncoor[connec[0], 0]; y = ncoor[connec[0], 1]
    gp = 1 / mp.sqrt(3)
    xsi = np.array([-gp, gp, gp, -gp])
    eta = np.array([-gp, -gp, gp, gp])

    eps_ip_1 = np.zeros(nb_inc_total)
    sig_ip_1 = np.zeros(nb_inc_total)
    phi_max  = np.zeros(nb_inc_total)
    energy_pos = np.zeros(nb_inc_total)

    # Historique & irréversibilité
    H_ip = np.zeros((4, 1))
    phi_node_old = np.zeros((4, 1))

    for it in range(nb_inc_total):
        u2 = u2_history[it]
        # DOF: [u1,v1,u2,v2,u3,v3,u4,v4]
        U = np.array([0.0, 0.0, u2, 0.0, 0.0, 0.0, 0.0, 0.0])

        Psi0_ip = np.zeros((4, 1))

        # --- Énergies positives par PG ---
        for ip in range(4):
            B_ip, _ = Strain(xsi[ip], eta[ip], x, y)
            eps = B_ip @ U  # [exx, eyy, gxy]
            if mode == "oed":
                eps_plus, eps_minus = oed_split_eps(eps, method=oed_method)
                Psi0_ip[ip, 0] = 0.5 * (eps_plus @ (C @ eps_plus))
            else:  # classic
                Psi0_ip[ip, 0] = classic_energy_positive(eps)

            if ip == 0:
                eps_ip_1[it] = eps[0]
                energy_pos[it] = Psi0_ip[ip, 0]

        # --- Historique H (monotone) ---
        H_ip = np.maximum(H_ip, Psi0_ip)

        # --- Champ de phase (H) + projection nodale d'irréversibilité ---
        K_phi = K_phi_integ(x, y, H_ip)
        R_phi = Resi_integ(x, y, H_ip)
        phi_node = np.linalg.solve(K_phi, R_phi.T)
        phi_node = np.maximum(phi_node, phi_node_old)
        phi_node = np.clip(phi_node, 0.0, 1.0)
        phi_node_old = phi_node.copy()
        phi_max[it] = np.max(phi_node)

        # --- φ aux PG ---
        def phi_at_ip(xsi_v, eta_v, phi_n):
            N = np.array([[0.25*(1-xsi_v)*(1-eta_v),
                           0.25*(1+xsi_v)*(1-eta_v),
                           0.25*(1+xsi_v)*(1+eta_v),
                           0.25*(1-xsi_v)*(1+eta_v)]])
            return (N @ phi_n)[0, 0]
        phi_ipp = np.zeros((4, 1))
        for ip in range(4):
            phi_ipp[ip, 0] = phi_at_ip(xsi[ip], eta[ip], phi_node)

        # --- Contraintes ---
        for ip in range(4):
            B_ip, _ = Strain(xsi[ip], eta[ip], x, y)
            eps = B_ip @ U
            if mode == "oed":
                eps_plus, eps_minus = oed_split_eps(eps, method=oed_method)
                g_d = (1.0 - phi_ipp[ip, 0])**2 + k_res
                stress = g_d * (C @ eps_plus) + (C @ eps_minus)
            else:
                g_d = (1.0 - phi_ipp[ip, 0])**2 + k_res
                stress = g_d * (C @ eps)   # dégradation de la contrainte totale
            if ip == 0:
                sig_ip_1[it] = stress[0]

    return {
        "eps_ip_1": eps_ip_1,
        "sig_ip_1": sig_ip_1,
        "phi_max":  phi_max,
        "energy_pos": energy_pos,
        "u2_history": u2_history.copy(),
        "mode": mode,
        "oed_method": oed_method
    }

# ------------- Lancement: OED(projector) vs OED(eig) vs Classic -------------
res_oed_proj = run_simulation(mode="oed", oed_method="projector")
res_oed_eig  = run_simulation(mode="oed", oed_method="eig")
res_classic  = run_simulation(mode="classic", oed_method="projector")  # oed_method ignoré ici

# ------------- Export des données en .dat -------------
increments = np.arange(nb_inc_total)

# Export 1: Chargement
with open('chargement.dat', 'w') as f:
    f.write("# Increment  u2_imposed\n")
    for i in range(nb_inc_total):
        f.write(f"{increments[i]:4d}  {res_oed_proj['u2_history'][i]:12.6f}\n")

# Export 2: σ-ε
with open('sigma_epsilon.dat', 'w') as f:
    f.write("# epsilon11_OED_proj  sigma11_OED_proj  epsilon11_OED_eig  sigma11_OED_eig  epsilon11_classic  sigma11_classic\n")
    for i in range(nb_inc_total):
        f.write(f"{res_oed_proj['eps_ip_1'][i]:12.6f}  {res_oed_proj['sig_ip_1'][i]:12.6f}  ")
        f.write(f"{res_oed_eig['eps_ip_1'][i]:12.6f}  {res_oed_eig['sig_ip_1'][i]:12.6f}  ")
        f.write(f"{res_classic['eps_ip_1'][i]:12.6f}  {res_classic['sig_ip_1'][i]:12.6f}\n")

# Export 3: φ max
with open('phi_max.dat', 'w') as f:
    f.write("# Increment  phi_max_OED_proj  phi_max_OED_eig  phi_max_classic\n")
    for i in range(nb_inc_total):
        f.write(f"{increments[i]:4d}  {res_oed_proj['phi_max'][i]:12.6f}  ")
        f.write(f"{res_oed_eig['phi_max'][i]:12.6f}  {res_classic['phi_max'][i]:12.6f}\n")

# Export 4: Ψ+
with open('energy_pos.dat', 'w') as f:
    f.write("# Increment  Psi_pos_OED_proj  Psi_pos_OED_eig  Psi_pos_classic\n")
    for i in range(nb_inc_total):
        f.write(f"{increments[i]:4d}  {res_oed_proj['energy_pos'][i]:12.6f}  ")
        f.write(f"{res_oed_eig['energy_pos'][i]:12.6f}  {res_classic['energy_pos'][i]:12.6f}\n")

print("Fichiers .dat exportés: chargement.dat, sigma_epsilon.dat, phi_max.dat, energy_pos.dat")

# ------------- Tracés comparatifs -------------
fig, axs = plt.subplots(2, 2, figsize=(14, 10))

# 1) σ-ε (PointGauss1)
ax = axs[0, 1]
ax.plot(res_oed_proj["eps_ip_1"], res_oed_proj["sig_ip_1"], '-o', ms=3, label='OED projector')
ax.plot(res_oed_eig["eps_ip_1"],  res_oed_eig["sig_ip_1"],  '--', ms=3, label='OED eig')
ax.plot(res_classic["eps_ip_1"],  res_classic["sig_ip_1"],  '-^', ms=3, label='Classic')
ax.set_xlabel("ε11"); ax.set_ylabel("σ11")
ax.set_title("σ–ε (PointGauss1) – comparaison")
ax.grid(alpha=0.3); ax.legend()
add_zero_axes(ax)

# 2) Chargement imposé u2
ax = axs[0, 0]
ax.plot(increments, res_oed_proj["u2_history"], '-', lw=2, color='orange')
ax.set_xlabel("Incrément")
ax.set_ylabel("u2 imposé")
ax.set_title("Chargement triangulaire imposé")
ax.grid(alpha=0.3)
ax.axhline(0, color='k', ls='--', lw=1, alpha=0.6)
ax.axvline(n1, color='b', ls=':', lw=2, alpha=0.8, label='Inversions')
ax.axvline(n1+n2, color='b', ls=':', lw=2, alpha=0.8)
ax.axvline(n1+n2+n3, color='b', ls=':', lw=2, alpha=0.8)
ax.legend()

# 3) φ max
axs[1, 0].plot(res_oed_proj["phi_max"], '-o',ms=3, label='OED projector')
axs[1, 0].plot(res_oed_eig["phi_max"],  '--',ms=3, label='OED eig')
axs[1, 0].plot(res_classic["phi_max"],  '-^',ms=3, label='Classic')
axs[1, 0].set_xlabel("Incrément"); axs[1, 0].set_ylabel("φ max")
axs[1, 0].set_title("Évolution du dommage")
axs[1, 0].grid(alpha=0.3); axs[1, 0].legend()

# 4) Ψ+
axs[1, 1].plot(res_oed_proj["energy_pos"], '-o',ms=3, label='OED projector')
axs[1, 1].plot(res_oed_eig["energy_pos"],  '--',ms=3, label='OED eig')
axs[1, 1].plot(res_classic["energy_pos"],  '-^',ms=3, alpha=0.5, label='Classic')
axs[1, 1].set_xlabel("Incrément"); axs[1, 1].set_ylabel("Ψ⁺ (PG1)")
axs[1, 1].set_title("Énergie élastique positive")
axs[1, 1].grid(alpha=0.3); axs[1, 1].legend()

plt.tight_layout()
plt.show()

print("\nRésumé :")
print(f" OED projector -> φ_max(final) = {res_oed_proj['phi_max'][-1]:.6e}, σ11_max = {np.max(res_oed_proj['sig_ip_1']):.6e}")
print(f" OED eig       -> φ_max(final) = {res_oed_eig['phi_max'][-1]:.6e}, σ11_max = {np.max(res_oed_eig['sig_ip_1']):.6e}")
print(f" Classic       -> φ_max(final) = {res_classic['phi_max'][-1]:.6e}, σ11_max = {np.max(res_classic['sig_ip_1']):.6e}")
