"""Shared likelihoods: cosmic chronometers, DESI DR2 BAO, Pantheon+ and DES-SN5YR (Dovekie)."""
import numpy as np
from scipy.integrate import cumulative_trapezoid

C_KMS = 299792.458
ZG = np.linspace(0, 2.34, 23401)


def comoving(E):
    with np.errstate(all="ignore"):
        e = E(ZG)
    if not np.all(np.isfinite(e)) or np.any(e <= 0):
        return None
    return cumulative_trapezoid(1 / e, ZG, initial=0)


def transverse(chi, ok):
    if ok > 1e-8:
        return np.sinh(np.sqrt(ok) * chi) / np.sqrt(ok)
    if ok < -1e-8:
        return np.sin(np.sqrt(-ok) * chi) / np.sqrt(-ok)
    return chi


# ---- cosmic chronometers (Moresco et al., BC03, total errors) ----
CC = np.array([
    [0.1791, 75.0, 4.0], [0.1993, 75.0, 5.0], [0.3519, 83.0, 14.0], [0.3802, 83.0, 13.5],
    [0.4004, 77.0, 10.2], [0.4247, 87.1, 11.2], [0.4497, 92.8, 12.9], [0.4783, 80.9, 9.0],
    [0.5929, 104.0, 13.0], [0.6797, 92.0, 8.0], [0.7812, 105.0, 12.0], [0.8754, 125.0, 17.0],
    [1.037, 154.0, 20.0], [1.363, 160.0, 33.6], [1.965, 186.5, 50.4]])


def chi2_cc(E, H0):
    with np.errstate(all="ignore"):
        m = H0 * E(CC[:, 0])
    if not np.all(np.isfinite(m)):
        return np.inf
    return float(np.sum(((CC[:, 1] - m) / CC[:, 2]) ** 2))


# ---- DESI DR2 BAO (Abdul Karim et al. 2025), as tabulated in arXiv:2604.11106 Table 2 ----
BGS = (0.295, 7.942, 0.075)
ANISO = [(0.510, 13.588, 0.167, 21.863, 0.425, -0.459), (0.706, 17.351, 0.177, 19.455, 0.330, -0.404),
         (0.934, 21.576, 0.152, 17.641, 0.193, -0.416), (1.321, 27.601, 0.318, 14.176, 0.221, -0.434),
         (1.484, 30.512, 0.760, 12.817, 0.516, -0.500), (2.330, 38.988, 0.531, 8.632, 0.101, -0.431)]
BAO_BLOCKS = [(z, np.array([dm, dh]), np.linalg.inv(np.array([[sdm**2, r * sdm * sdh], [r * sdm * sdh, sdh**2]])))
              for z, dm, sdm, dh, sdh, r in ANISO]


def chi2_bao(E, ok, H0, rd):
    chi = comoving(E)
    if chi is None:
        return np.inf
    dm_all = transverse(chi, ok) * C_KMS / H0 / rd
    with np.errstate(all="ignore"):
        dh = lambda z: C_KMS / (H0 * E(z)) / rd
    zb, dv, sdv = BGS
    dmb = np.interp(zb, ZG, dm_all)
    c2 = ((dv - (zb * dmb**2 * dh(zb)) ** (1 / 3)) / sdv) ** 2
    for z, d, icov in BAO_BLOCKS:
        r = d - np.array([np.interp(z, ZG, dm_all), dh(z)])
        c2 += r @ icov @ r
    return float(c2)


# ---- supernovae ----
class SN:
    def __init__(self, which="pantheon"):
        if which == "pantheon":
            tab = np.genfromtxt("data/Pantheon+SH0ES.dat", names=True, dtype=None, encoding=None)
            with open("data/Pantheon+SH0ES_STAT+SYS.cov") as f:
                n = int(f.readline())
                cov = np.loadtxt(f).reshape(n, n)
            m = tab["zHD"] > 0.01
            self.z, self.zhel, self.mu = tab["zHD"][m], tab["zHEL"][m], tab["m_b_corr"][m]
            self.icov = np.linalg.inv(cov[np.ix_(m, m)])
            self.name = f"Pantheon+ ({m.sum()} SNe)"
        else:
            folder = "data/des_repo/4_DISTANCES_COVMAT/"
            rows, header = [], None
            for line in open(folder + "DES-Dovekie_HD.csv"):
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.replace(",", " ").split()
                if parts[0] in ("VARNAMES:", "SN:"):
                    parts = parts[1:]
                if header is None:
                    header = parts
                else:
                    rows.append(parts)
            col = {k: i for i, k in enumerate(header)}
            get = lambda k: np.array([float(r[col[k]]) for r in rows])
            self.z, self.zhel, self.mu = get("zHD"), get("zHEL"), get("MU")
            d = np.load(folder + "STAT+SYS.npz")
            n = int(d["nsn"][0])
            icov = np.zeros((n, n))
            icov[np.triu_indices(n)] = d["cov"]
            il = np.tril_indices(n, -1)
            icov[il] = icov.T[il]
            self.icov = icov
            self.name = f"DES-SN5YR Dovekie ({n} SNe)"
        self.ones = np.ones(len(self.z))
        self.S = self.ones @ self.icov @ self.ones
        if self.z.max() > ZG.max():
            raise ValueError("supernova redshifts beyond the distance grid")

    def chi2(self, E, ok):
        chi = comoving(E)
        if chi is None:
            return np.inf
        with np.errstate(all="ignore"):
            mu = 5 * np.log10((1 + self.zhel) * np.interp(self.z, ZG, transverse(chi, ok)))
        if not np.all(np.isfinite(mu)):
            return np.inf
        r = self.mu - mu
        B = self.ones @ self.icov @ r
        return float(r @ self.icov @ r - B**2 / self.S)
