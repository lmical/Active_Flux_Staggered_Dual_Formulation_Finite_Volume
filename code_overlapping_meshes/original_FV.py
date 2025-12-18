#!/usr/bin/env python3
"""
Second-order finite-volume (MUSCL + minmod, loop-based) with Rusanov flux and SSPRK3.
Solves 1D Euler equations (example: Sod shock tube).
Boundary conditions: independent left/right, with outflow, wall, inflow options.
"""

import numpy as np
import matplotlib.pyplot as plt

# --------------------------
# Utility: slope limiter
# --------------------------
def minmod(a, b):
    """Scalar minmod limiter."""
    if a * b <= 0.0:
        return 0.0
    else:
        return np.sign(a) * min(abs(a), abs(b))

# --------------------------
# Euler equations
# --------------------------
gamma = 1.4

def primitive_from_cons(U):
    """Convert conserved to primitive. U: (3,) array."""
    rho = U[0]
    u = U[1] / rho
    E = U[2]
    p = (gamma - 1.0) * (E - 0.5 * rho * u**2)
    return rho, u, p

def cons_from_prim(rho, u, p):
    """Convert primitive to conserved variables."""
    E = p / (gamma - 1.0) + 0.5 * rho * u**2
    return np.array([rho, rho * u, E])

def flux_euler(U):
    """Physical flux for Euler equations. U: (3,) array."""
    rho, u, p = primitive_from_cons(U)
    return np.array([rho * u, rho * u**2 + p, u * (U[2] + p)])

def max_wave_speed(U):
    """Max wave speed |u|+c. U: (3,) array."""
    rho, u, p = primitive_from_cons(U)
    c = np.sqrt(max(gamma * p / rho, 1e-12))
    return abs(u) + c

# --------------------------
# Boundary conditions
# --------------------------
def apply_bc(U, left_bc=("outflow",), right_bc=("outflow",)):
    """
    Apply boundary conditions with 2 ghost cells.
    Each bc is a tuple: ("outflow",) or ("wall",) or ("inflow", rho,u,p).
    """
    # --- Left boundary ---
    bc = left_bc[0]
    if bc == "outflow":
        U[:, 0] = U[:, 2]
        U[:, 1] = U[:, 2]
    elif bc == "wall":
        U[:, 0] = U[:, 2]; U[:, 1] = U[:, 2]
        U[1, 0] = -U[1, 2]; U[1, 1] = -U[1, 2]  # reverse momentum
    elif bc == "inflow":
        rho, u, p = left_bc[1], left_bc[2], left_bc[3]
        UL = cons_from_prim(rho, u, p)
        U[:, 0] = UL; U[:, 1] = UL
    else:
        raise ValueError(f"Unknown left BC: {bc}")

    # --- Right boundary ---
    bc = right_bc[0]
    if bc == "outflow":
        U[:, -1] = U[:, -3]
        U[:, -2] = U[:, -3]
    elif bc == "wall":
        U[:, -1] = U[:, -3]; U[:, -2] = U[:, -3]
        U[1, -1] = -U[1, -3]; U[1, -2] = -U[1, -3]
    elif bc == "inflow":
        rho, u, p = right_bc[1], right_bc[2], right_bc[3]
        UR = cons_from_prim(rho, u, p)
        U[:, -1] = UR; U[:, -2] = UR
    else:
        raise ValueError(f"Unknown right BC: {bc}")

# --------------------------
# MUSCL reconstruction & Rusanov flux
# --------------------------
def reconstruct_MUSCL(U):
    """
    MUSCL reconstruction with minmod limiter (loop-based).
    Input: U (3, N) with ghosts.
    Output: UL, UR at interfaces (3, N-3).
    """
    nv, N = U.shape
    UL = np.zeros((nv, N-3))
    UR = np.zeros((nv, N-3))

    # loop over interior cells (1..N-2)
    for i in range(1, N-1):
        for k in range(nv):
            # compute slope for cell i
            dl = U[k, i] - U[k, i-1]
            dr = U[k, i+1] - U[k, i]
            slope = minmod(dl, dr)

            # left state at i+1/2 comes from cell i
            if i < N-2:
                UL[k, i-1] = U[k, i] + 0.5 * slope

            # right state at i-1/2 comes from cell i
            if i > 1:
                UR[k, i-2] = U[k, i] - 0.5 * slope

    return UL, UR

def rusanov_flux(UL, UR):
    """Compute Rusanov flux at interfaces (loop-based)."""
    nv, M = UL.shape
    F = np.zeros((nv, M))
    for j in range(M):
        ULj = UL[:, j]
        URj = UR[:, j]
        FL = flux_euler(ULj)
        FR = flux_euler(URj)
        a = max(max_wave_speed(ULj), max_wave_speed(URj))
        F[:, j] = 0.5 * (FL + FR) - 0.5 * a * (URj - ULj)
    return F

# --------------------------
# RHS and time integrator
# --------------------------
def rhs(U, dx, left_bc, right_bc):
    nv, N = U.shape
    apply_bc(U, left_bc=left_bc, right_bc=right_bc)
    UL, UR = reconstruct_MUSCL(U)
    F = rusanov_flux(UL, UR)
    dUdt = np.zeros_like(U)

    # update interior cells (exclude 2 ghost cells on each side)
    for i in range(2, N-2):
        dUdt[:, i] = -(F[:, i-1] - F[:, i-2]) / dx
    return dUdt

def ssp_rk3_step(U, dt, dx, left_bc, right_bc):
    U0 = U.copy()
    k1 = rhs(U, dx, left_bc, right_bc)
    U1 = U0 + dt * k1
    apply_bc(U1, left_bc, right_bc)

    k2 = rhs(U1, dx, left_bc, right_bc)
    U2 = 0.75 * U0 + 0.25 * (U1 + dt * k2)
    apply_bc(U2, left_bc, right_bc)

    k3 = rhs(U2, dx, left_bc, right_bc)
    Unew = (1.0/3.0) * U0 + (2.0/3.0) * (U2 + dt * k3)
    apply_bc(Unew, left_bc, right_bc)

    return Unew

# --------------------------
# Initial condition: Sod shock tube
# --------------------------
def sod_initial(x):
    U = np.zeros((3, x.size))
    for i in range(x.size):
        if np.abs(x[i]-0.5)<1e-12:
            U[0, i]=0.5*(1.0+0.125)
            U[1, i]=0.0
            U[2, i]=0.5*(1.0+0.1)/(gamma-1)
        elif (x[i]<0.5):
            rho = 1.0
            u = 0.0
            p = 1.0
            U[:, i] = cons_from_prim(rho, u, p)
        else:
            rho = 0.125
            u = 0.0
            p = 0.1
            U[:, i] = cons_from_prim(rho, u, p)
    return U

# --------------------------
# Driver
# --------------------------
def run_sod(nx=200, cfl=0.475, t_final=0.2,
            left_bc=("outflow",), right_bc=("outflow",)):
    xL, xR = 0.0, 1.0
    dx = (xR - xL) / nx
    xc = np.linspace(xL + 0.5*dx, xR - 0.5*dx, nx)

    U_phys = sod_initial(xc)
    U = np.zeros((3, nx + 4))
    U[:, 2:-2] = U_phys
    apply_bc(U, left_bc, right_bc)

    t = 0.0
    indt=0
    while t < t_final:
        indt=indt+1
        # compute dt from CFL
        amax = 0.0
        for i in range(2, nx+2):
            amax = max(amax, max_wave_speed(U[:, i]))
        dt = cfl * dx / amax
        print("Iteration", indt,"from time",f"{t:.14f}","with dt",f"{dt:.14f}")
        if t + dt > t_final:
            dt = t_final - t
        U = ssp_rk3_step(U, dt, dx, left_bc, right_bc)
        t += dt

    return U[:, 2:-2], xc

# --------------------------
# Main
# --------------------------
if __name__ == "__main__":
    # Example: Sod with outflow BCs
    U, x = run_sod(nx=200, cfl=0.475, t_final=0.2,
                   left_bc=("outflow",), right_bc=("outflow",))

    rho, u, p = [], [], []
    for i in range(U.shape[1]):
        ri, ui, pi = primitive_from_cons(U[:, i])
        rho.append(ri); u.append(ui); p.append(pi)

    with open("solution", "w") as f:
        f.write("# x rho rhou E\n")
        nx=len(x)
        for i in range(nx):
            f.write(f"{x[i]:.14f} {U[0,i]:.14f} {U[1,i]:.14f} {U[2,i]:.14f}\n")



    plt.figure(figsize=(4,6))
    plt.subplot(311); plt.plot(x, rho, "-k"); plt.ylabel("rho")
    plt.subplot(312); plt.plot(x, u, "-k");   plt.ylabel("u")
    plt.subplot(313); plt.plot(x, p, "-k");   plt.ylabel("p"); plt.xlabel("x")
    plt.tight_layout()
    plt.show()
