#!/usr/bin/env python
"""Primordial abundances z=1000->10, v2026 Grackle (Peebles C-factor
recombination + CMB Compton + CMB H-/H2+ dissociation, the latter now in the
correct CODE units = *time_units).  Three panels:
  (1) all species fractions, with the RECFAST x_e overlay;
  (2) H- and H2+ with CMB dissociation ON vs OFF (shows the terms are now applied
      -- ON pins them to the floor where the CMB destroys them);
  (3) temperature: T_gas now Compton-locks to T_CMB at high z and decouples ~z150
      (no spurious below-CMB cooling once H2 is at its true ~1e-6 level)."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

C = dict(z=0,Tg=1,Tc=2,xHII=3,xH2=4,xHM=5,xH2p=6,xe=7,xHD=8,xDII=9)
d1 = np.loadtxt("/tmp/sp_diss1.dat")   # CMB dissociation ON  (the corrected run)
d0 = np.loadtxt("/tmp/sp_diss0.dat")   # CMB dissociation OFF
z, zp1 = d1[:,0], 1+d1[:,0]
rf = np.loadtxt("/Users/tabel/Projects/cicass/vbc_transfer/recfast/xeTrecfast.out", skiprows=1)
m = (rf[:,0] >= 10) & (rf[:,0] <= 1000)

fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(7.4, 11.4), sharex=True,
                                 gridspec_kw=dict(height_ratios=[2.0,1.5,1.0]))

# (1) species
a1.loglog(zp1, d1[:,C['xe']], "C3", lw=2.2, label=r"$x_e$ (Grackle v2026)")
a1.loglog(1+rf[m,0], rf[m,1], "k--", lw=1.5, label=r"$x_e$ (RECFAST)")
a1.loglog(zp1, d1[:,C['xH2']], "C0", lw=2.0, label=r"$x_{\rm H_2}$")
a1.loglog(zp1, d1[:,C['xHD']], "C2", lw=2.0, label=r"$x_{\rm HD}$")
a1.loglog(zp1, d1[:,C['xDII']],"C5", lw=1.4, ls="-.", label=r"$x_{\rm D^+}$")
a1.set_ylim(1e-13, 3e-1); a1.set_ylabel(r"n(X)/n$_{\rm H}$")
a1.legend(ncol=2, fontsize=9, loc="lower left"); a1.grid(which="both", alpha=0.15)
a1.set_title("Primordial chemistry z=1000$\\to$10 — Grackle v2026 "
             "(C-factor recomb. + CMB dissociation)", fontsize=10.5)
a1.annotate(r"$x_{\rm H_2}\!\to\!3\times10^{-6}$ (canonical)", (12, 3.1e-6),
            (30, 2e-5), fontsize=8.5, arrowprops=dict(arrowstyle="->", lw=0.8))

# (2) H- and H2+ : dissociation ON vs OFF
a2.loglog(zp1, d1[:,C['xHM']],  "C1", lw=2.2, label=r"$x_{\rm H^-}$  CMB diss. ON")
a2.loglog(1+d0[:,0], d0[:,C['xHM']], "C1", lw=1.5, ls=":", label=r"$x_{\rm H^-}$  OFF")
a2.loglog(zp1, d1[:,C['xH2p']], "C4", lw=2.2, label=r"$x_{\rm H_2^+}$  CMB diss. ON")
a2.loglog(1+d0[:,0], d0[:,C['xH2p']],"C4", lw=1.5, ls=":", label=r"$x_{\rm H_2^+}$  OFF")
a2.set_ylim(1e-21, 1e-1); a2.set_ylabel(r"n(X)/n$_{\rm H}$")
a2.legend(ncol=2, fontsize=8.5, loc="lower right"); a2.grid(which="both", alpha=0.15)
a2.text(0.02,0.06,"CMB photo-destruction now applied in code units\n"
        "(*time_units) → H$^-$, H$_2^+$ pinned to the floor at high z",
        transform=a2.transAxes, fontsize=8.2,
        bbox=dict(boxstyle="round", fc="#eef7ee", ec="0.6"))

# (3) temperature
a3.loglog(zp1, d1[:,C['Tg']], "C3", lw=2.0, label=r"$T_{\rm gas}$ (Grackle, diss. ON)")
a3.loglog(zp1, d0[:,C['Tg']], "C0", lw=1.4, ls=":", label=r"$T_{\rm gas}$ (diss. OFF)")
a3.loglog(zp1, d1[:,C['Tc']], "0.4", lw=1.5, ls="--", label=r"$T_{\rm CMB}$")
a3.loglog(1+rf[m,0], rf[m,2], "k:", lw=1.3, label=r"$T_{\rm gas}$ (RECFAST)")
a3.set_ylabel("T [K]"); a3.set_xlabel("1 + z")
a3.legend(fontsize=8.3, loc="upper left"); a3.grid(which="both", alpha=0.15)
a3.set_xlim(zp1.max()*1.05, zp1.min()*0.95)

fig.tight_layout()
out = "/Users/tabel/Research/codes/grackle/src/example/species_v2026.png"
fig.savefig(out, dpi=140); print("wrote", out)
