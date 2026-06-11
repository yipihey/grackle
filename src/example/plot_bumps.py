#!/usr/bin/env python
"""The two H2-formation bumps (Galli & Palla 1998): as the universe expands the
CMB stops photo-dissociating H2+ (first bump, z~400-600) and then stops photo-
detaching H- (second bump, z~80-150), so the two H2-formation channels switch on
in sequence. Computed from the H2-formation rate through each channel with the
CMB survival factor (using OUR fork's GP98 CMB rates k27/k28) on the recfast/CICASS
recombination history n_e(z), T(z)."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

# recfast/CICASS history
R = np.loadtxt("/Users/tabel/Projects/cicass/vbc_transfer/recfast/xeTrecfast.out", skiprows=1)
rz, rxe, rT = R[:,0], R[:,1], R[:,2]
o = np.argsort(rz); rz, rxe, rT = rz[o], rxe[o], rT[o]

z = np.linspace(800, 10, 1600)
xe = np.interp(z, rz, rxe); Tg = np.interp(z, rz, rT)
Tr = 2.73*(1+z)                                   # CMB temperature
h, Ob, XH = 0.71, 0.046, 0.76
nH = Ob*1.8788e-29*h*h*XH/1.6726e-24*(1+z)**3     # physical n_H [cm^-3]
ne = xe*nH; nHII = xe*nH

# rate coefficients (cm^3/s); Abel+1997 / Galli&Palla 1998
k7  = 3.0e-16*(Tg/300.)**0.95*np.exp(-Tg/9320.)   # H + e -> H- + g
k8  = 1.3e-9*np.ones_like(Tg)                     # H- + H -> H2 + e
k9  = 1.85e-23*Tg**1.8                            # H+ + H -> H2+ + g (radiative assoc)
k10 = 6.0e-10*np.ones_like(Tg)                    # H2+ + H -> H2 + H+
# OUR fork's CMB photo-rates (s^-1), GP98
k27 = 1.1e-1*Tr**2.13*np.exp(-8823./Tr)           # H- + g_CMB -> H + e
k28 = 1.63e7*np.exp(-32400./Tr)                   # H2+ + g_CMB -> H + H+ (LTE)

# H2 formation rate per H nucleus through each channel = (formation) x (CMB-survival)
surv_HM  = k8*nH /(k8*nH  + k27)                  # fraction of H- that makes H2
surv_H2p = k10*nH/(k10*nH + k28)                  # fraction of H2+ that makes H2
R_HM  = k7*ne *surv_HM  * nH / nH                 # = k7 ne surv_HM   [s^-1 per H]
R_H2p = k9*nHII*surv_H2p* nH / nH                 # = k9 nHII surv_H2p
# (per-H rates: R_HM = k7*ne*surv_HM, R_H2p = k9*nHII*surv_H2p)
R_HM  = k7*ne*surv_HM
R_H2p = k9*nHII*surv_H2p

# cumulative H2 fraction: integrate dx(H2)/dt = R over cosmic time (EdS, z>10)
H0 = h*100*1e5/3.0857e24; Om=0.27
t = (2./3.)/H0/np.sqrt(Om)*(1+z)**-1.5            # cosmic time [s], increasing as z drops
xH2 = np.concatenate([[0], np.cumsum(0.5*(R_HM+R_H2p)[1:]+0.5*(R_HM+R_H2p)[:-1])*np.diff(t)])
xH2_HM  = np.concatenate([[0], np.cumsum(0.5*(R_HM)[1:]+0.5*(R_HM)[:-1])*np.diff(t)])
xH2_H2p = np.concatenate([[0], np.cumsum(0.5*(R_H2p)[1:]+0.5*(R_H2p)[:-1])*np.diff(t)])

fig,(ax0,ax1)=plt.subplots(2,1,figsize=(7.2,8.6),sharex=True)
# --- the two bumps in the formation rate ---
ax0.semilogy(1+z, R_H2p, "b-", lw=2.2, label=r"H$_2^+$ channel  (k9·n$_{\rm H^+}$·surv$_{\rm H_2^+}$)")
ax0.semilogy(1+z, R_HM,  "r-", lw=2.2, label=r"H$^-$ channel  (k7·n$_e$·surv$_{\rm H^-}$)")
ax0.semilogy(1+z, R_H2p+R_HM, "k-", lw=1.0, alpha=0.6, label="total")
ax0.set_ylabel(r"H$_2$ formation rate per H  [s$^{-1}$]")
ax0.set_ylim(1e-22,3e-17); ax0.legend(frameon=False,fontsize=9,loc="lower center")
ax0.set_title("Two H$_2$-formation bumps as the CMB releases H$_2^+$ then H$^-$ (our CMB rates)")
# mark the two release redshifts (where survival crosses 0.5)
for surv, c, lab in [(surv_H2p,"b","H$_2^+$ release"),(surv_HM,"r","H$^-$ release")]:
    i=np.argmin(np.abs(surv-0.5)); ax0.axvline(1+z[i],color=c,ls=":",lw=1.2)
    ax0.text(1+z[i], 1.5e-17, f"z≈{z[i]:.0f}", color=c, ha="center", fontsize=8)

# --- cumulative H2 (two-step rise) ---
ax1.semilogy(1+z, np.maximum(xH2_H2p,1e-12), "b--", lw=1.6, label="via H$_2^+$")
ax1.semilogy(1+z, np.maximum(xH2_HM,1e-12),  "r--", lw=1.6, label="via H$^-$")
ax1.semilogy(1+z, np.maximum(xH2,1e-12),     "k-",  lw=2.4, label="total H$_2$")
ax1.set_xlabel(r"$1+z$"); ax1.set_ylabel(r"$n_{\rm H_2}/n_{\rm H}$ (cumulative)")
ax1.set_ylim(1e-9,1e-5); ax1.legend(frameon=False,fontsize=9,loc="upper right")
for ax in (ax0,ax1):
    ax.set_xlim((1+z).max(),(1+z).min()); ax.grid(True,which="major",alpha=0.25)
fig.tight_layout(); fig.savefig("/tmp/h2_two_bumps.png",dpi=140)
print("wrote /tmp/h2_two_bumps.png")
ih2p=np.argmax(R_H2p); ihm=np.argmax(R_HM)
print(f"H2+ bump peak at z={z[ih2p]:.0f} (T_cmb={Tr[ih2p]:.0f}K); H- bump peak at z={z[ihm]:.0f} (T_cmb={Tr[ihm]:.0f}K)")
print(f"final n(H2)/n_H = {xH2[-1]:.2e}  (via H2+: {xH2_H2p[-1]:.2e}, via H-: {xH2_HM[-1]:.2e})")
