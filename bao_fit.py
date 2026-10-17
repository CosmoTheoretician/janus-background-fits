import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# DESI DR2 BAO (Abdul Karim et al. 2025, PRD 112, 083515), as tabulated in arXiv:2604.11106, Table 2
# tracer, z_eff, DM/rd, sigma, DH/rd, sigma, rho(DM, DH); BGS gives DV/rd only
BGS = (0.295, 7.942, 0.075)
ANISO = [
    ("LRG1", 0.510, 13.588, 0.167, 21.863, 0.425, -0.459),
    ("LRG2", 0.706, 17.351, 0.177, 19.455, 0.330, -0.404),
    ("LRG3+ELG1", 0.934, 21.576, 0.152, 17.641, 0.193, -0.416),
    ("ELG2", 1.321, 27.601, 0.318, 14.176, 0.221, -0.434),
    ("QSO", 1.484, 30.512, 0.760, 12.817, 0.516, -0.500),
    ("Lya", 2.330, 38.988, 0.531, 8.632, 0.101, -0.431),
]
blocks = []
for name, z, dm, sdm, dh, sdh, rho in ANISO:
    cov = np.array([[sdm**2, rho * sdm * sdh], [rho * sdm * sdh, sdh**2]])
    blocks.append((name, z, np.array([dm, dh]), np.linalg.inv(cov)))


def E_janus(z, q0):
    return (1 + z) * np.sqrt(1 + 2 * q0 * z)


def E_lcdm(z, om):
    return np.sqrt(om * (1 + z) ** 3 + 1 - om)


def E_milne(z):
    return 1 + z


def distances(z, E, ok):
    """DM and DH in units of c/H0, for curvature density ok (>0 open)."""
    chi = quad(lambda x: 1 / E(x), 0, z, limit=200)[0]
    if ok > 1e-8:
        dm = np.sinh(np.sqrt(ok) * chi) / np.sqrt(ok)
    else:
        dm = chi
    return dm, 1 / E(z)


def chi2(E, ok, alpha):
    """alpha = c / (H0 rd), the only scale BAO can constrain."""
    c2 = 0.0
    zb, dv, sdv = BGS
    dm, dh = distances(zb, E, ok)
    c2 += ((dv - alpha * (zb * dm**2 * dh) ** (1 / 3)) / sdv) ** 2
    for name, z, d, icov in blocks:
        dm, dh = distances(z, E, ok)
        r = d - alpha * np.array([dm, dh])
        c2 += r @ icov @ r
    return c2


QMIN = -0.999 / (2 * 2.330)
models = {
    "Janus (q0 free)": (lambda p: chi2(lambda z: E_janus(z, p[1]), 1 - 2 * p[1], p[0]), [30, -0.05], [(10, 60), (QMIN, 2)]),
    "Janus (q0 = -0.087)": (lambda p: chi2(lambda z: E_janus(z, -0.087), 1 + 0.174, p[0]), [30], [(10, 60)]),
    "Flat LCDM (Om free)": (lambda p: chi2(lambda z: E_lcdm(z, p[1]), 0.0, p[0]), [30, 0.3], [(10, 60), (0.01, 0.99)]),
    "Milne": (lambda p: chi2(E_milne, 1.0, p[0]), [30], [(10, 60)]),
}
n = 13
lines = ["DESI DR2 BAO: 13 measurements, 0.295 < z < 2.33, r_d free (only c/(H0 rd) fitted)", "",
         f"{'model':24s} {'chi2':>8s} {'dof':>4s} {'chi2/dof':>9s} {'AIC':>8s} {'BIC':>8s}  best fit [c/(H0 rd), shape]"]
best = {}
for name, (f, p0, bounds) in models.items():
    res = minimize(f, p0, bounds=bounds, method="L-BFGS-B")
    k = len(p0)
    best[name] = res
    lines.append(f"{name:24s} {res.fun:8.2f} {n-k:4d} {res.fun/(n-k):9.2f} {res.fun+2*k:8.2f} {res.fun+k*np.log(n):8.2f}  {np.round(res.x, 4)}")

q_grid = np.linspace(QMIN + 1e-3, 0.5, 120)
prof = np.array([minimize(lambda p: chi2(lambda z: E_janus(z, q), 1 - 2 * q, p[0]), [30],
                          bounds=[(10, 60)], method="L-BFGS-B").fun for q in q_grid])
qb = q_grid[np.argmin(prof)]
inside = q_grid[prof <= prof.min() + 1]
lines += ["", f"Janus q0 profile: best {qb:.3f}, 68% range {inside.min():.3f} to {inside.max():.3f}",
          f"chi2(q0 = -0.087) - chi2(best Janus) = {np.interp(-0.087, q_grid, prof) - prof.min():.1f}",
          f"chi2(best Janus) - chi2(best flat LCDM) = {best['Janus (q0 free)'].fun - best['Flat LCDM (Om free)'].fun:.1f}"]
text = "\n".join(lines)
print(text)
open("results/bao_fit.txt", "w").write(text + "\n")

# plot: data divided by best-fit flat LCDM prediction
a_l, om = best["Flat LCDM (Om free)"].x
ref = lambda z: distances(z, lambda x: E_lcdm(x, om), 0.0)
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
zz = np.linspace(0.2, 2.4, 60)
curves = {"Janus (q0 free)": "C3", "Janus (q0 = -0.087)": "C1", "Milne": "C2"}
for i, lab in enumerate(["D_M / r_d", "D_H / r_d"]):
    ax = axes[i]
    for name, z, d, icov in blocks:
        r = ref(z)[i] * a_l
        err = np.sqrt(np.linalg.inv(icov)[i, i])
        ax.errorbar(z, d[i] / r, err / r, fmt="o", color="k")
    for name, col in curves.items():
        p = best[name].x
        q = p[1] if len(p) > 1 else (-0.087 if "0.087" in name else None)
        if name == "Milne":
            E, ok = E_milne, 1.0
        else:
            E, ok = (lambda x, q=q: E_janus(x, q)), 1 - 2 * q
        ax.plot(zz, [p[0] * distances(x, E, ok)[i] / (a_l * ref(x)[i]) for x in zz], color=col, label=name)
    ax.axhline(1, color="C0", label="flat LCDM (best fit)")
    ax.set_xlabel("z"); ax.set_ylabel(f"{lab}, relative to best flat LCDM"); ax.legend()
fig.savefig("results/bao_fit.png", dpi=120, bbox_inches="tight")
fig.savefig("/mnt/c/Users/roche/Desktop/janus_bao_fit.png", dpi=120, bbox_inches="tight")
print("\nresults saved to results/bao_fit.txt and results/bao_fit.png (copy on your desktop)")
