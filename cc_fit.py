import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Moresco et al. cosmic chronometers (BC03), total uncertainties, km/s/Mpc
data = np.array([
    [0.1791, 75.0, 4.0], [0.1993, 75.0, 5.0], [0.3519, 83.0, 14.0], [0.3802, 83.0, 13.5],
    [0.4004, 77.0, 10.2], [0.4247, 87.1, 11.2], [0.4497, 92.8, 12.9], [0.4783, 80.9, 9.0],
    [0.5929, 104.0, 13.0], [0.6797, 92.0, 8.0], [0.7812, 105.0, 12.0], [0.8754, 125.0, 17.0],
    [1.037, 154.0, 20.0], [1.363, 160.0, 33.6], [1.965, 186.5, 50.4]])
z, H, sig = data.T
zmax = z.max()


def H_janus(z, H0, q0):
    arg = 1 + 2 * q0 * z
    return np.where(arg > 0, H0 * (1 + z) * np.sqrt(np.clip(arg, 0, None)), np.nan)


def H_lcdm(z, H0, om):
    return H0 * np.sqrt(om * (1 + z) ** 3 + 1 - om)


def H_milne(z, H0):
    return H0 * (1 + z)


def chi2(model, p):
    m = model(z, *p)
    return np.inf if np.any(~np.isfinite(m)) else float(np.sum(((H - m) / sig) ** 2))


fits = {}
specs = {
    "Janus (H0, q0 free)": (H_janus, [70, -0.1], [(40, 100), (-0.999 / (2 * zmax), 2)]),
    "Janus (q0 = -0.087)": (lambda z, H0: H_janus(z, H0, -0.087), [70], [(40, 100)]),
    "Flat LCDM (H0, Om free)": (H_lcdm, [70, 0.3], [(40, 100), (0.0, 1.0)]),
    "Milne (H0 free)": (H_milne, [70], [(40, 100)]),
}
n = len(z)
lines = [f"{n} cosmic chronometer points, 0.18 < z < {zmax}", "",
         f"{'model':28s} {'chi2':>7s} {'dof':>4s} {'chi2/dof':>9s} {'AIC':>7s} {'BIC':>7s}  best fit"]
for name, (model, p0, bounds) in specs.items():
    res = minimize(lambda p: chi2(model, p), p0, bounds=bounds, method="L-BFGS-B")
    k = len(p0)
    c2 = res.fun
    fits[name] = (model, res.x)
    lines.append(f"{name:28s} {c2:7.2f} {n-k:4d} {c2/(n-k):9.2f} {c2+2*k:7.2f} {c2+k*np.log(n):7.2f}  {np.round(res.x, 3)}")

q_grid = np.linspace(-0.25, 0.6, 171)
prof = np.array([minimize(lambda p: chi2(H_janus, [p[0], q]), [70], bounds=[(40, 100)], method="L-BFGS-B").fun
                 for q in q_grid])
best = q_grid[np.argmin(prof)]
inside = q_grid[prof <= prof.min() + 1]
lines += ["", f"Janus q0 profile: best {best:.3f}, 68% range {inside.min():.3f} to {inside.max():.3f}",
          f"chi2 at q0 = -0.087 minus best: {np.interp(-0.087, q_grid, prof) - prof.min():.2f}"]
text = "\n".join(lines)
print(text)
open("results/cc_fit.txt", "w").write(text + "\n")

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
zz = np.linspace(0, 2.1, 300)
axes[0].errorbar(z, H, sig, fmt="o", color="k", label="cosmic chronometers")
for (name, (model, p)), c in zip(fits.items(), ["C3", "C1", "C0", "C2"]):
    axes[0].plot(zz, model(zz, *p), color=c, label=name)
axes[0].set_xlabel("z"); axes[0].set_ylabel("H(z) [km/s/Mpc]"); axes[0].legend()
axes[1].plot(q_grid, prof - prof.min(), color="C3")
axes[1].axvline(-0.087, ls="--", color="C1", label="Petit & d'Agostini q0 = -0.087")
axes[1].axhline(1, ls=":", color="gray")
axes[1].set_xlabel("q0"); axes[1].set_ylabel("Δχ² (Janus)"); axes[1].legend()
fig.savefig("results/cc_fit.png", dpi=120, bbox_inches="tight")
fig.savefig("/mnt/c/Users/roche/Desktop/janus_cc_fit.png", dpi=120, bbox_inches="tight")
print("\nresults saved to results/cc_fit.txt and results/cc_fit.png (copy on your desktop)")
