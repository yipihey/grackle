#!/usr/bin/env python
"""Primordial species fractions z=1000 -> 10: the v2026 MINIMAL reduced network
(advect only HII, H2I, HDI -- He neutral, e=p, H-/H2+/D+ equilibrium, D & HI
reconstructed) overlaid on the FULL 12-species Grackle network (all species
advected).  Same CMB recombination + CMB dissociation physics.  The minimal
model carries 3 fields instead of 12 yet tracks the full network to <0.1%."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

# cols: z T_gas T_cmb xHII xH2 xHM xH2p xe xHD xDII
F = np.loadtxt("/tmp/species_full.dat")
M = np.loadtxt("/tmp/species_min.dat")
zF, zM = 1+F[:,0], 1+M[:,0]
def col(A, i): return np.clip(A[:,i], 1e-30, None)

fig, (ax, axr) = plt.subplots(2, 1, figsize=(7.2, 8.2), sharex=True,
                              gridspec_kw=dict(height_ratios=[3, 1]))

species = [(3,"$x_{\\rm HII}=x_e$","C3"), (4,"$x_{\\rm H_2}$","C0"),
           (8,"$x_{\\rm HD}$","C2"), (5,"$x_{\\rm H^-}$","C1"),
           (6,"$x_{\\rm H_2^+}$","C4"), (9,"$x_{\\rm D^+}$","C5")]
# the species columns in the file: index 3=xHII,4=xH2,5=xHM,6=xH2p,7=xe,8=xHD,9=xDII
idx = {"HII":3,"H2":4,"HM":5,"H2p":6,"xe":7,"HD":8,"DII":9}
plot = [("HII","$x_{\\rm HII}$","C3"),("H2","$x_{\\rm H_2}$","C0"),
        ("HD","$x_{\\rm HD}$","C2"),("HM","$x_{\\rm H^-}$","C1"),
        ("H2p","$x_{\\rm H_2^+}$","C4"),("DII","$x_{\\rm D^+}$","C5")]
for key,lab,c in plot:
    j = idx[key]
    ax.loglog(zF, col(F,j), color=c, lw=2.0, label=lab+" (full)")
    ax.loglog(zM[::12], col(M,j)[::12], color=c, ls="none", marker="o", ms=3.5,
              mfc="none", label=lab+" (minimal)")
ax.set_ylim(1e-14, 3e-1); ax.set_ylabel("abundance fraction n(X)/n$_{\\rm H}$")
ax.legend(ncol=3, fontsize=7.6, loc="lower left", framealpha=0.9)
ax.set_title("Minimal (3 advected: HII,H2I,HDI) vs full 12-species network, "
             "z=1000$\\to$10", fontsize=10.5)
ax.grid(which="both", alpha=0.15)

# relative error panel for the three ADVECTED species
for key,lab,c in [("HII","HII","C3"),("H2","H2","C0"),("HD","HD","C2")]:
    j = idx[key]
    rel = np.abs(col(M,j)-col(F,j))/np.maximum(col(F,j),1e-30)
    axr.loglog(zM, np.maximum(rel,1e-7), color=c, lw=1.6, label=lab)
axr.axhline(1e-2, color="0.6", ls=":", lw=1); axr.text(900, 1.2e-2, "1%", fontsize=8, color="0.4")
axr.set_ylim(1e-6, 1e-1); axr.set_ylabel("|minimal-full|/full"); axr.set_xlabel("1 + z")
axr.legend(fontsize=8, loc="upper left", ncol=3); axr.grid(which="both", alpha=0.15)
axr.set_xlim(zF.max()*1.05, zF.min()*0.95)

fig.tight_layout()
out = "/Users/tabel/Research/codes/grackle/src/example/minimal_vs_full.png"
fig.savefig(out, dpi=140); print("wrote", out)
for key,j in [("HII",3),("H2",4),("HD",8)]:
    m = zM < 300
    rel = np.max(np.abs(col(M,j)[m]-col(F,j)[m])/col(F,j)[m])
    print(f"  max |minimal-full|/full for {key} (z<300): {rel:.2e}")
