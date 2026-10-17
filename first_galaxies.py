import numpy as np
from scipy.integrate import quad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GYR = 977.8  # (1/H0) in Gyr for H0 in km/s/Mpc
GALAXIES = [("GN-z11", 10.6), ("JADES-GS-z13-0", 13.2), ("JADES-GS-z14-0", 14.2), ("MoM-z14", 14.4)]
OLDEST_STARS = (12.5, 13.5)  # rough range for the oldest globular clusters / halo stars, Gyr

E_janus = lambda q: (lambda z: (1 + z) * np.sqrt(1 + 2 * q * z))
E_lcdm = lambda om: (lambda z: np.sqrt(om * (1 + z) ** 3 + 1 - om))
# H0 values for Janus taken from the combined-fit profiles (q0 = 0 -> 61.2, -0.087 -> 63.7, +0.085 -> 59.0)
MODELS = [
    ("Flat LCDM (combined fit)", E_lcdm(0.304), 68.1, None),
    ("Janus q0 = -0.087 (Petit), H0 = 63.7", E_janus(-0.087), 63.7, 1 / (2 * 0.087)),
    ("Janus q0 = -0.048 (growth best), H0 ~ 62.5", E_janus(-0.048), 62.5, 1 / (2 * 0.048)),
    ("Janus q0 = +0.085 (combined best), H0 = 59.0", E_janus(0.085), 59.0, None),
]


def age(E, H0, z, zmax):
    """time from z_max (or the big bang) to redshift z, in Gyr."""
    top = np.inf if zmax is None else zmax * (1 - 1e-9)
    if zmax is not None and z >= zmax:
        return 0.0
    return GYR / H0 * quad(lambda x: 1 / ((1 + x) * E(x)), z, top, limit=400)[0]


lines = [f"Cosmic ages (Gyr). Oldest stars: about {OLDEST_STARS[0]}-{OLDEST_STARS[1]} Gyr.", ""]
hdr = f"{'model':46s} {'age today':>9s} {'z_start':>8s}" + "".join(f" {n+' z='+str(z):>22s}" for n, z in GALAXIES)
lines.append(hdr)
fig, ax = plt.subplots(figsize=(8, 5))
zz = np.logspace(-2, np.log10(30), 400)
for name, E, H0, zmax in MODELS:
    t0 = age(E, H0, 0, zmax)
    row = f"{name:46s} {t0:9.2f} {('inf' if zmax is None else f'{zmax:.2f}'):>8s}"
    for gname, zg in GALAXIES:
        if zmax is not None and zg >= zmax:
            row += f" {'before matter era':>22s}"
        else:
            row += f" {age(E, H0, zg, zmax):22.3f}"
    lines.append(row)
    zp = zz if zmax is None else zz[zz < zmax]
    ax.plot(zp, [age(E, H0, z, zmax) for z in zp], label=name)
    if zmax is not None:
        ax.plot([zmax], [0], "x", color=ax.lines[-1].get_color())
lines += ["", "For Janus with q0 < 0, 'age' counts from the start of the dust-era solution at z_start = 1/(2|q0|),",
          "where the scale factor is at its minimum; the model places an earlier radiative era before that,",
          "with different physics, so galaxies at z > z_start fall outside the solution that was fitted."]
text = "\n".join(lines)
print(text)
open("results/first_galaxies.txt", "w").write(text + "\n")

for gname, zg in GALAXIES:
    ax.axvline(zg, color="gray", ls=":", lw=0.8)
ax.text(GALAXIES[-1][1] * 1.05, 10, "JWST galaxies", color="gray", fontsize=8)
ax.axhspan(*OLDEST_STARS, color="gold", alpha=0.25, label="oldest stars")
ax.set_xscale("log"); ax.set_xlabel("z"); ax.set_ylabel("time since start [Gyr]"); ax.legend(fontsize=7)
fig.savefig("results/first_galaxies.png", dpi=120, bbox_inches="tight")
fig.savefig("/mnt/c/Users/roche/Desktop/janus_first_galaxies.png", dpi=120, bbox_inches="tight")
print("\nresults saved to results/first_galaxies.txt and results/first_galaxies.png")
