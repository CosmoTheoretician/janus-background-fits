import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.optimize import minimize_scalar

FOLDER = "data/des_repo/4_DISTANCES_COVMAT/"

# Hubble diagram: skip '#' comment lines, then a header row (possibly 'VARNAMES:'), then data rows (possibly 'SN:')
rows, header = [], None
for line in open(FOLDER + "DES-Dovekie_HD.csv"):
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    parts = [p for p in line.replace(",", " ").split()]
    if parts[0] in ("VARNAMES:", "SN:"):
        parts = parts[1:]
    if header is None:
        header = parts
    else:
        rows.append(parts)
col = {name: i for i, name in enumerate(header)}
get = lambda name: np.array([float(r[col[name]]) for r in rows])
z, zhel, mu = get("zHD"), get("zHEL"), get("MU")

# inverse covariance (stat+sys), stored as the upper triangle, as in DES-Dovekie-SN_Likelihood.py
d = np.load(FOLDER + "STAT+SYS.npz")
n = int(d["nsn"][0])
icov = np.zeros((n, n))
icov[np.triu_indices(n)] = d["cov"]
il = np.tril_indices(n, -1)
icov[il] = icov.T[il]
assert n == len(z), f"covariance size {n} != {len(z)} supernovae"
keep = z > 0
icov = icov[np.ix_(keep, keep)]  # all kept in DES's own likelihood; kept here for safety
z, zhel, mu = z[keep], zhel[keep], mu[keep]
ones = np.ones(len(z))
S = ones @ icov @ ones
n = len(z)
zg = np.linspace(0, z.max() * 1.01, 20001)


def mu_shape(E, ok):
    with np.errstate(over="ignore"):
        chi = cumulative_trapezoid(1 / E(zg), zg, initial=0)
        dm = np.sinh(np.sqrt(ok) * chi) / np.sqrt(ok) if ok > 1e-8 else chi
    return 5 * np.log10((1 + zhel) * np.interp(z, zg, dm))


def chi2(E, ok):
    m = mu_shape(E, ok)
    if not np.all(np.isfinite(m)):
        return np.inf
    r = mu - m
    B = ones @ icov @ r
    return r @ icov @ r - B**2 / S


E_janus = lambda q: (lambda x: (1 + x) * np.sqrt(np.clip(1 + 2 * q * x, 1e-12, None)))
E_lcdm = lambda om: (lambda x: np.sqrt(om * (1 + x) ** 3 + 1 - om))
qmin = -0.999 / (2 * z.max())
rj = minimize_scalar(lambda q: chi2(E_janus(q), 1 - 2 * q), bounds=(qmin, 1.0), method="bounded")
rl = minimize_scalar(lambda om: chi2(E_lcdm(om), 0.0), bounds=(0.01, 0.99), method="bounded")
c_petit = chi2(E_janus(-0.087), 1.174)
c_milne = chi2(lambda x: 1 + x, 1.0)
q_grid = np.linspace(max(qmin, -0.3), 0.3, 121)
prof = np.array([chi2(E_janus(q), 1 - 2 * q) for q in q_grid])
inside = q_grid[prof <= prof.min() + 1]

text = "\n".join([
    f"DES-SN5YR (Dovekie recalibration): {n} supernovae ({z.min():.3f} < z < {z.max():.2f}),",
    "full stat+sys covariance, absolute magnitude marginalised", "",
    f"{'model':26s} {'chi2':>9s} {'dof':>5s} {'chi2/dof':>9s}  best fit",
    f"{'Janus (q0 free)':26s} {rj.fun:9.2f} {n-2:5d} {rj.fun/(n-2):9.3f}  q0 = {rj.x:.3f}",
    f"{'Janus (q0 = -0.087)':26s} {c_petit:9.2f} {n-1:5d} {c_petit/(n-1):9.3f}",
    f"{'Flat LCDM (Om free)':26s} {rl.fun:9.2f} {n-2:5d} {rl.fun/(n-2):9.3f}  Om = {rl.x:.3f}",
    f"{'Milne':26s} {c_milne:9.2f} {n-1:5d} {c_milne/(n-1):9.3f}", "",
    f"Janus q0 68% range: {inside.min():.3f} to {inside.max():.3f}",
    f"chi2(best Janus) - chi2(best flat LCDM) = {rj.fun - rl.fun:.1f}"])
print(text)
open("results/des_fit.txt", "w").write(text + "\n")
