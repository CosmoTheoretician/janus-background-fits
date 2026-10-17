import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize_scalar
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

tab = np.genfromtxt("data/Pantheon+SH0ES.dat", names=True, dtype=None, encoding=None)
with open("data/Pantheon+SH0ES_STAT+SYS.cov") as f:
    n_all = int(f.readline())
    cov_all = np.loadtxt(f).reshape(n_all, n_all)
mask = tab["zHD"] > 0.01
z, zhel, mb = tab["zHD"][mask], tab["zHEL"][mask], tab["m_b_corr"][mask]
cov = cov_all[np.ix_(mask, mask)]
icov = np.linalg.inv(cov)
ones = np.ones(len(z))
S = ones @ icov @ ones
n = len(z)

zg = np.linspace(0, z.max() * 1.01, 20001)


def mu_shape(E, ok):
    chi = cumulative_trapezoid(1 / E(zg), zg, initial=0)
    dm = np.sinh(np.sqrt(ok) * chi) / np.sqrt(ok) if ok > 1e-8 else chi
    return 5 * np.log10((1 + zhel) * np.interp(z, zg, dm))


def chi2(E, ok):
    """chi2 with the magnitude offset (M and H0) marginalised analytically."""
    if not np.all(np.isfinite(E(zg))):
        return np.inf
    r = mb - mu_shape(E, ok)
    B = ones @ icov @ r
    return r @ icov @ r - B**2 / S


E_janus = lambda q: (lambda x: (1 + x) * np.sqrt(np.clip(1 + 2 * q * x, 1e-12, None)))
E_lcdm = lambda om: (lambda x: np.sqrt(om * (1 + x) ** 3 + 1 - om))
qmin = -0.999 / (2 * z.max())

rj = minimize_scalar(lambda q: chi2(E_janus(q), 1 - 2 * q), bounds=(qmin, 1.0), method="bounded")
rl = minimize_scalar(lambda om: chi2(E_lcdm(om), 0.0), bounds=(0.01, 0.99), method="bounded")
c_petit = chi2(E_janus(-0.087), 1 + 0.174)
c_milne = chi2(lambda x: 1 + x, 1.0)

lines = [f"Pantheon+ (z > 0.01): {n} supernovae, full STAT+SYS covariance, absolute magnitude marginalised", "",
         f"{'model':26s} {'chi2':>9s} {'dof':>5s} {'chi2/dof':>9s}  best fit",
         f"{'Janus (q0 free)':26s} {rj.fun:9.2f} {n-2:5d} {rj.fun/(n-2):9.3f}  q0 = {rj.x:.3f}",
         f"{'Janus (q0 = -0.087)':26s} {c_petit:9.2f} {n-1:5d} {c_petit/(n-1):9.3f}",
         f"{'Flat LCDM (Om free)':26s} {rl.fun:9.2f} {n-2:5d} {rl.fun/(n-2):9.3f}  Om = {rl.x:.3f}",
         f"{'Milne':26s} {c_milne:9.2f} {n-1:5d} {c_milne/(n-1):9.3f}"]
q_grid = np.linspace(max(qmin, -0.3), 0.3, 121)
prof = np.array([chi2(E_janus(q), 1 - 2 * q) for q in q_grid])
inside = q_grid[prof <= prof.min() + 1]
lines += ["", f"Janus q0 68% range: {inside.min():.3f} to {inside.max():.3f}",
          f"chi2(best Janus) - chi2(best flat LCDM) = {rj.fun - rl.fun:.1f}"]
text = "\n".join(lines)
print(text)
open("results/sn_fit.txt", "w").write(text + "\n")

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
ref = mu_shape(E_lcdm(rl.x), 0.0)
off_ref = (ones @ icov @ (mb - ref)) / S
order = np.argsort(z)
bins = np.logspace(np.log10(0.01), np.log10(z.max()), 20)
idx = np.digitize(z, bins)
zb = [z[idx == i].mean() for i in range(1, len(bins)) if np.any(idx == i)]
rb = [np.mean((mb - ref - off_ref)[idx == i]) for i in range(1, len(bins)) if np.any(idx == i)]
axes[0].scatter(z, mb - ref - off_ref, s=2, color="0.7")
axes[0].plot(zb, rb, "ko", label="binned data")
for (name, E, ok, col) in [("Janus (q0 free)", E_janus(rj.x), 1 - 2 * rj.x, "C3"),
                          ("Janus (q0 = -0.087)", E_janus(-0.087), 1.174, "C1"),
                          ("Milne", lambda x: 1 + x, 1.0, "C2")]:
    m = mu_shape(E, ok)
    off = (ones @ icov @ (mb - m)) / S
    axes[0].plot(z[order], (m + off - ref - off_ref)[order], color=col, label=name)
axes[0].axhline(0, color="C0", label="flat LCDM")
axes[0].set_xscale("log"); axes[0].set_ylim(-0.4, 0.4)
axes[0].set_xlabel("z"); axes[0].set_ylabel("distance modulus relative to best flat LCDM"); axes[0].legend()
axes[1].plot(q_grid, prof - prof.min(), color="C3")
axes[1].axvline(-0.087, ls="--", color="C1", label="Petit & d'Agostini q0 = -0.087")
axes[1].set_xlabel("q0"); axes[1].set_ylabel("Δχ² (Janus, Pantheon+)"); axes[1].legend()
fig.savefig("results/sn_fit.png", dpi=120, bbox_inches="tight")
fig.savefig("/mnt/c/Users/roche/Desktop/janus_sn_fit.png", dpi=120, bbox_inches="tight")
print("\nresults saved to results/sn_fit.txt and results/sn_fit.png (copy on your desktop)")
