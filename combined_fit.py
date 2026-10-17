import sys
import numpy as np
from scipy.optimize import minimize
import emcee
import corner
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from likelihoods import SN, chi2_cc, chi2_bao, CC

which = sys.argv[1] if len(sys.argv) > 1 else "pantheon"
sn = SN(which)
N_DATA = len(CC) + 13 + len(sn.z)

MODELS = {
    "Janus": dict(
        names=["H0", "q0", "rd"], p0=[65, -0.03, 150],
        bounds=[(50, 90), (-0.999 / (2 * 2.34), 0.6), (100, 200)],
        bg=lambda p: ((lambda z, q=p[1]: (1 + z) * np.sqrt(np.clip(1 + 2 * q * z, 1e-15, None))), 1 - 2 * p[1])),
    "Flat LCDM": dict(
        names=["H0", "Om", "rd"], p0=[68, 0.3, 147],
        bounds=[(50, 90), (0.05, 0.7), (100, 200)],
        bg=lambda p: ((lambda z, om=p[1]: np.sqrt(om * (1 + z) ** 3 + 1 - om)), 0.0)),
    "oLCDM (curvature + Lambda)": dict(
        names=["H0", "Om", "OL", "rd"], p0=[68, 0.3, 0.7, 147],
        bounds=[(50, 90), (0.05, 0.7), (0.0, 1.2), (100, 200)],
        bg=lambda p: ((lambda z, om=p[1], ol=p[2]: np.sqrt(om * (1 + z) ** 3 + (1 - om - ol) * (1 + z) ** 2 + ol)),
                      1 - p[1] - p[2])),
}


def parts(model, p):
    E, ok = model["bg"](p)
    H0, rd = p[0], p[-1]
    return chi2_cc(E, H0), chi2_bao(E, ok, H0, rd), sn.chi2(E, ok)


def total(model, p):
    if any(not (lo <= x <= hi) for x, (lo, hi) in zip(p, model["bounds"])):
        return np.inf
    return sum(parts(model, p))


lines = [f"Combined fit: 15 cosmic chronometers + 13 DESI DR2 BAO + {sn.name}; r_d free, SN magnitude marginalised", "",
         f"{'model':28s} {'chi2_CC':>8s} {'chi2_BAO':>9s} {'chi2_SN':>9s} {'total':>9s} {'k':>2s} {'AIC':>9s} {'BIC':>9s}"]
best = {}
for name, m in MODELS.items():
    res = minimize(lambda p: total(m, p), m["p0"], method="Nelder-Mead", options=dict(maxiter=20000, xatol=1e-6, fatol=1e-6))
    best[name] = res.x
    cc, bao, s = parts(m, res.x)
    k = len(res.x)
    lines.append(f"{name:28s} {cc:8.2f} {bao:9.2f} {s:9.2f} {res.fun:9.2f} {k:2d} {res.fun+2*k:9.2f} {res.fun+k*np.log(N_DATA):9.2f}")
    lines.append(f"{'':28s} best fit: " + ", ".join(f"{a} = {b:.4g}" for a, b in zip(m["names"], res.x)))

lines += ["", "MCMC posteriors (median and 68% interval):"]
for name in ["Janus", "Flat LCDM", "oLCDM (curvature + Lambda)"]:
    m = MODELS[name]
    ndim, nwalk = len(m["names"]), 24
    start = best[name] + 1e-3 * np.random.default_rng(1).standard_normal((nwalk, ndim)) * np.abs(best[name])
    sampler = emcee.EnsembleSampler(nwalk, ndim, lambda p: -0.5 * total(m, p))
    sampler.run_mcmc(start, 2500, progress=True)
    chain = sampler.get_chain(discard=700, thin=5, flat=True)
    desc = []
    for i, pn in enumerate(m["names"]):
        lo, med, hi = np.percentile(chain[:, i], [16, 50, 84])
        desc.append(f"{pn} = {med:.4g} (+{hi-med:.2g} / -{med-lo:.2g})")
    lines.append(f"  {name}: " + ", ".join(desc))
    fig = corner.corner(chain, labels=m["names"], show_titles=True)
    tag = name.split()[0].lower()
    fig.savefig(f"results/combined_corner_{tag}_{which}.png", dpi=100)
    plt.close(fig)

dj = total(MODELS["Janus"], best["Janus"]) - total(MODELS["Flat LCDM"], best["Flat LCDM"])
lines += ["", f"chi2(Janus) - chi2(flat LCDM) = {dj:.1f}  (same number of free parameters)"]
text = "\n".join(lines)
print("\n" + text)
open(f"results/combined_fit_{which}.txt", "w").write(text + "\n")
