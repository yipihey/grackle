#!/usr/bin/env python
"""Galli & Palla (1998)-style figure: primordial species fractions + thermal
history vs redshift, from the Grackle cosmology one-zone (cosmo_evolve.c,
z=200 -> 10 (post-recombination)) using our fork's CMB photo-rates."""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.loadtxt("/tmp/cosmo_species.dat")
z, Tg, Tcmb, xHII, xH2, xHM, xH2p, xe, xHD, xDII = d.T
zp1 = 1.0 + z

fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(7.0, 8.4), sharex=True,
                               gridspec_kw=dict(height_ratios=[1, 1.5]))

# --- thermal history ---
ax0.loglog(zp1, Tcmb, "k--", lw=1.6, label=r"$T_{\rm CMB}=2.73(1+z)$")
ax0.loglog(zp1, Tg, "r-", lw=2.0, label=r"$T_{\rm gas}$ (Grackle: Compton + adiabatic)")
ax0.set_ylabel("temperature  [K]")
ax0.set_ylim(1, 5e3)
ax0.legend(frameon=False, fontsize=9, loc="lower right")
ax0.set_title("Primordial chemistry vs redshift (Grackle network; recfast/CICASS recombination ICs)")

# --- species fractions (Galli & Palla 1998 style) ---
sp = [(xe,   "e$^-$ , H$^+$",     "k-",  2.0),
      (xH2,  "H$_2$",            "b-",  2.0),
      (xHD,  "HD",               "g-",  1.6),
      (xHM,  "H$^-$",            "m-",  1.4),
      (xH2p, "H$_2^+$",          "c-",  1.4),
      (xDII, "D$^+$",            "orange", 1.4)]
for y, lab, st, lw in sp:
    if isinstance(st, str) and st in ("orange",):
        ax1.loglog(zp1, np.maximum(y, 1e-30), color=st, lw=lw, label=lab)
    else:
        ax1.loglog(zp1, np.maximum(y, 1e-30), st, lw=lw, label=lab)
ax1.set_xlabel(r"$1+z$")
ax1.set_ylabel(r"fractional abundance  $n_X/n_{\rm H}$")
ax1.set_ylim(1e-14, 1.0)
ax1.set_xlim(zp1.max(), zp1.min())          # high-z (early) on the left
ax0.set_xlim(zp1.max(), zp1.min())
ax1.legend(frameon=False, fontsize=9, ncol=2, loc="lower left")
for ax in (ax0, ax1):
    ax.grid(True, which="major", alpha=0.25)
    for zz in (1000, 100, 10):
        ax.axvline(1+zz, color="grey", lw=0.5, alpha=0.4)

fig.tight_layout()
out = "/tmp/gp98_grackle.png"
fig.savefig(out, dpi=140)
print("wrote", out)
# quick numeric summary
for zz in (1000, 300, 100, 30, 10):
    i = np.argmin(np.abs(z - zz))
    print(f"z={z[i]:6.1f}  Tgas={Tg[i]:7.1f}  xe={xe[i]:.2e}  xH2={xH2[i]:.2e}  xHD={xHD[i]:.2e}")
