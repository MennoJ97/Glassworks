"""An independent reference for lamellar gratings: the Fourier modal method (RCWA).

It shares nothing with the Rust solver: no mesh, no absorbing layer, and fields expanded in plane
waves. E_z uses plain Fourier factorization, and H_z uses Li's inverse rule. Run with numpy:
    python3 rcwa.py
"""
import numpy as np


def layer_fourier(n_ridge, n_gap, x0, x1, period, size):
    """Toeplitz matrices of eps and of 1/eps, for a ridge from x0 to x1, gap elsewhere."""
    def coefficient(inside, outside, m):
        if m == 0:
            return outside + (inside - outside) * (x1 - x0) / period
        g = 2 * np.pi * m / period
        return (inside - outside) * (np.exp(-1j * g * x0) - np.exp(-1j * g * x1)) / (1j * g * period)
    def toeplitz(inside, outside):
        return np.array([[coefficient(inside, outside, i - j) for j in range(size)] for i in range(size)])
    return toeplitz(n_ridge**2, n_gap**2), toeplitz(1 / n_ridge**2, 1 / n_gap**2)


def grating(wavelength, angle_deg, period, n_top, n_bottom, n_ridge, n_gap, x0, x1, depth, polarization, orders=40):
    """Each propagating order's share of the incident power, as {m: share} on each side."""
    k0 = 2 * np.pi / wavelength
    ms = np.arange(-orders, orders + 1)
    size = len(ms)
    kx = k0 * n_top * np.sin(np.radians(angle_deg)) + 2 * np.pi * ms / period
    Kx = np.diag(kx.astype(complex))
    e, ei = layer_fourier(n_ridge, n_gap, x0, x1, period, size)
    if polarization == 'Ez':
        A = Kx @ Kx - k0**2 * e
        weight = np.eye(size, dtype=complex)
        top_p, bottom_p = 1.0, 1.0
    else:
        A = np.linalg.inv(ei) @ (Kx @ np.linalg.inv(e) @ Kx - k0**2 * np.eye(size))
        weight = ei
        top_p, bottom_p = 1 / n_top**2, 1 / n_bottom**2
    values, W = np.linalg.eig(A)
    # u'' = A u, so u = W exp(+-gamma y), with gamma on the branch that decays or goes outgoing.
    gamma = np.sqrt(values.astype(complex))
    gamma = np.where(gamma.real < 0, -gamma, gamma)
    ky = lambda n: np.sqrt((k0 * n)**2 - kx**2 + 0j)
    k_top, k_bottom = ky(n_top), ky(n_bottom)
    k_top = np.where(k_top.imag < 0, -k_top, k_top)
    k_bottom = np.where(k_bottom.imag < 0, -k_bottom, k_bottom)
    # Unknowns: r, then c+ and c- in the layer (with exponentials scaled), then t. Layer from y=0 to -depth.
    V = weight @ W @ np.diag(gamma)
    phase = np.diag(np.exp(-gamma * depth))
    Z = np.zeros((size, size), dtype=complex)
    I = np.eye(size, dtype=complex)
    # u = inc + r, at the top: W (c+ + phase c-) ; p du/dy: i ky p (r - inc) = V (c+ - phase c-)
    # At the bottom: W (phase c+ + c-) = t ; V (phase c+ - c-) = -i ky p t
    M = np.block([
        [-I, W, W @ phase, Z],
        [-1j * top_p * np.diag(k_top), V, -V @ phase, Z],
        [Z, W @ phase, W, -I],
        [Z, V @ phase, -V, 1j * bottom_p * np.diag(k_bottom)],
    ])
    incident = np.zeros(size, dtype=complex)
    incident[orders] = 1
    rhs = np.concatenate([incident, -1j * top_p * k_top * incident, np.zeros(2 * size, dtype=complex)])
    x = np.linalg.solve(M, rhs)
    r, t = x[:size], x[3 * size:]
    power = top_p * k_top[orders].real
    shares = lambda amp, k, p: {int(m): float(abs(a)**2 * p * kk.real / power) for m, a, kk in zip(ms, amp, k) if abs(kk.imag) < 1e-12 and kk.real > 0}
    return shares(r, k_top, top_p), shares(t, k_bottom, bottom_p)


if __name__ == '__main__':
    cases = [
        (1.0, 0.0, 1.5, 1.0, 1.5, 1.5, 1.0, 0.375, 1.125, 0.5),
        (1.0, 20.0, 1.5, 1.0, 1.5, 1.5, 1.0, 0.375, 1.125, 0.5),
        (1.0, 20.0, 1.5, 1.0, 1.5, 2.5, 1.0, 0.375, 0.975, 0.5),
    ]
    for case in cases:
        for pol in ['Ez', 'Hz']:
            print(case, pol)
            for orders in [40, 80, 160]:
                r, t = grating(*case, pol, orders=orders)
                total = sum(r.values()) + sum(t.values())
                print(f'  orders {orders}: total {total:.6f}')
                if orders == 160:
                    print('   R', {m: round(v, 5) for m, v in r.items()})
                    print('   T', {m: round(v, 5) for m, v in t.items()})
