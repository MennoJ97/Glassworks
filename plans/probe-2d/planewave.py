"""An independent reference for the rod lattice of chapters 9 and 11: a plane-wave expansion.

It shares nothing with the Rust solver: no mesh, no elements, and fields expanded in plane waves
(Ho, Chan & Soukoulis 1990). A circle's Fourier coefficient is exact. E_z uses that matrix of eps
directly, and H_z uses its inverse, the inverse rule. Run with numpy:
    python3 planewave.py
"""
import numpy as np


def bessel_j1(x):
    """J_1 from its integral, (1/pi) times the integral of cos(t - x sin t) over 0 to pi."""
    t = (np.arange(400) + 0.5) * np.pi / 400
    return np.cos(t - np.multiply.outer(x, np.sin(t))).mean(axis=-1)


def eps_matrix(gx, gy, eps_rod, radius, cell, rods):
    """eps(G - G') for rods of `radius` at `rods`, in air, in a square cell of side `cell`."""
    dx, dy = np.subtract.outer(gx, gx), np.subtract.outer(gy, gy)
    g = np.hypot(dx, dy)
    fill = np.pi * radius**2 / cell**2
    gr = np.where(g == 0, 1.0, g * radius)
    # Many pairs share a |G - G'|, so evaluate J_1 once for each.
    values, index = np.unique(gr.round(12), return_inverse=True)
    shape = np.where(g == 0, 1.0, 2 * bessel_j1(values)[index].reshape(g.shape) / gr)
    where = sum(np.exp(-1j * (dx * x + dy * y)) for x, y in rods)
    return np.eye(len(gx)) + (eps_rod - 1) * fill * shape * where


def frequencies(k, polarization, eps_rod=12.0, radius=0.2, cell=1.0, rods=((0.0, 0.0),), waves=24, count=4):
    """The lowest `count` frequencies a / lambda at Bloch wavevector `k`, with lattice period 1."""
    m = np.arange(-waves, waves + 1)
    mx, my = [a.ravel() for a in np.meshgrid(m, m)]
    gx, gy = 2 * np.pi * mx / cell, 2 * np.pi * my / cell
    qx, qy = k[0] + gx, k[1] + gy
    eps = eps_matrix(gx, gy, eps_rod, radius, cell, rods)
    if polarization == 'Ez':
        # |k + G|^2 E = w^2 eps E: a generalized problem, made standard through eps's Cholesky factor.
        l = np.linalg.cholesky(eps)
        li = np.linalg.inv(l)
        op = li @ np.diag(qx**2 + qy**2) @ li.conj().T
    else:
        # (k + G).(k + G') [eps]^-1 H = w^2 H.
        op = (np.outer(qx, qx) + np.outer(qy, qy)) * np.linalg.inv(eps)
    w2 = np.linalg.eigvalsh((op + op.conj().T) / 2)[:count]
    return np.sqrt(np.maximum(w2, 0)) / (2 * np.pi)


if __name__ == '__main__':
    pi = np.pi
    print('Band edges, square lattice of eps 12 rods, r = 0.2a. Chapter 9, last row.')
    for waves in [12, 18, 24]:
        n = (2 * waves + 1)**2
        for pol in ['Ez', 'Hz']:
            m1 = frequencies((pi, pi), pol, waves=waves)[0]
            x2 = frequencies((pi, 0), pol, waves=waves)[1]
            print(f'  {n:5d} waves, {pol}: M1 {m1:.5f}, X2 {x2:.5f}')
    print('One rod removed from a 5 x 5 supercell, at Gamma, E_z. Chapter 11.')
    rods = [(i, j) for i in range(-2, 3) for j in range(-2, 3) if (i, j) != (0, 0)]
    for waves in [21, 27, 33]:
        n = (2 * waves + 1)**2
        f = frequencies((0, 0), 'Ez', cell=5.0, rods=rods, waves=waves, count=40)
        inside = [x for x in f if 0.29 < x < 0.41]
        print(f'  {n:5d} waves: modes in the gap {np.round(inside, 5)}')
