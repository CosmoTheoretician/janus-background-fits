import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

C_CM = 2.99792458e10
DECADE = 10 * 3.15576e7
MODELS = [
    ("Flat LCDM (H0 = 68.1, Om = 0.304)", 68.1, lambda z: np.sqrt(0.304 * (1 + z) ** 3 + 0.696), None),
    ("Janus q0 = -0.087 (Petit), H0 = 63.7", 63.7, lambda z: (1 + z) * np.sqrt(1 - 0.174 * z), 5.75),
    ("Janus q0 = +0.085 (combined best), H0 = 59.0", 59.0, lambda z: (1 + z) * np.sqrt(1 + 0.170 * z), None),
    ("Milne (coasting), H0 = 61", 61.0, lambda z: 1 + z, None),
]
ZS = [0.5, 1, 2, 3, 4, 5]


def drift(H0, E, z):  # cm/s per decade
    return C_CM * (H0 * 1e3 / 3.0857e22) * DECADE * (1 - E(z) / (1 + z))


lines = ["Redshift drift: velocity shift of a fixed source, cm/s per decade (positive = redshift increasing)", "",
         f"{'model':46s}" + "".join(f" {'z='+str(z):>8s}" for z in ZS)]
fig, ax = plt.subplots(figsize=(8, 5))
zz = np.linspace(0.01, 6, 400)
for name, H0, E, zmax in MODELS:
    lines.append(f"{name:46s}" + "".join(f" {'n/a' if zmax and z >= zmax else f'{drift(H0, E, z):+.2f}':>8s}" for z in ZS))
    zp = zz if zmax is None else zz[zz < zmax]
    ax.plot(zp, drift(H0, E, zp), label=name)
lines += ["", "ELT/ANDES forecast: a few cm/s precision over 10-20 years in the Lyman-alpha forest, z ~ 2-5."]
text = "\n".join(lines)
print(text)
open("results/redshift_drift.txt", "w").write(text + "\n")
ax.axhline(0, color="gray", lw=0.8)
ax.axvspan(2, 5, color="gold", alpha=0.2, label="ELT Lyα forest range")
ax.set_xlabel("z"); ax.set_ylabel("Δv per decade [cm/s]"); ax.legend(fontsize=8)
fig.savefig("results/redshift_drift.png", dpi=120, bbox_inches="tight")
fig.savefig("/mnt/c/Users/roche/Desktop/janus_redshift_drift.png", dpi=120, bbox_inches="tight")
