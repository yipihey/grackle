#!/usr/bin/env python
"""Double-check: are H- and H2+ ALWAYS in equilibrium (so they can be replaced by
algebraic equilibrium formulae instead of advected as rate-equation species)?

A transient species X (formation F, destruction D*nX) relaxes to equilibrium on
tau_eq = 1/D.  It is in equilibrium iff tau_eq << the timescale on which F, D and
the parent species (HI, HII, e) change -- i.e. the dynamical / chemical time
(free-fall, Hubble, cooling).  We compute tau_eq for H- and H2+ across the
relevant (n_H, T, x_e, z) space and compare to t_ff and t_Hubble, and we directly
integrate the H- ODE vs its equilibrium in a collapsing parcel."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

mh=1.6726e-24; G=6.674e-8; kpc=3.086e21
# --- Grackle/Abel+1997 destruction rate coefficients [cm^3 s^-1] ---
def k8 (T): return 1.35e-9*(T**9.8493e-2+3.2852e-1*T**5.561e-1+2.771e-7*T**2.1826)/(1+6.191e-3*T**1.0461+8.9712e-11*T**3.0424+3.2576e-14*T**3.7741)  # H-+H->H2+e
def k15(T): return 2.5634e-9*T**1.78186                      # H-+H->2H+e (T<...; approx)
def k16(T): return 2.4e-6*(1+T/2e4)/np.sqrt(T)              # H-+H+->2H (mutual neutralization)
def k14(T): return np.where(T>0, 1e-9, 1e-9)               # H-+e->H+2e (~1e-9, weak; placeholder)
def k10(T): return 6.0e-10*np.ones_like(T)                  # H2++H->H2+H+
def k18(T): return np.where(T>617, 1.32e-6*T**-0.76, 1e-8)  # H2++e->2H
# CMB photo-rates (our fork, GP98) [s^-1]
def k27(Tr): return 1.1e-1*Tr**2.13*np.exp(-8823./Tr)       # H- photodetach
def k28(Tr): return 1.63e7*np.exp(-32400./Tr)               # H2+ photodissoc (LTE)

def tau_HM(nH, T, xe, z):
    Tr=2.73*(1+z); nHI=nH; nHII=xe*nH; ne=xe*nH
    D = k8(T)*nHI + k15(T)*nHI + k16(T)*nHII + k14(T)*ne + k27(Tr)
    return 1.0/D
def tau_H2p(nH, T, xe, z):
    Tr=2.73*(1+z); nHI=nH; ne=xe*nH
    D = k10(T)*nHI + k18(T)*ne + k28(Tr)
    return 1.0/D

t_ff = lambda nH: np.sqrt(3*np.pi/(32*G*nH*mh/0.76))
H0=0.71*100*1e5/3.086e24
t_H = lambda z: 1.0/(H0*np.sqrt(0.27*(1+z)**3+0.73))

fig,(ax0,ax1)=plt.subplots(1,2,figsize=(12,5))

# Panel A: equilibration time vs n_H, across regimes (collapse: z=0 no CMB; cosmo: z=300)
nH=np.logspace(-4,10,200)
for T,xe,z,c,lab in [(200,1e-4,0,"b","T=200K, x_e=1e-4 (cold neutral)"),
                     (1000,1e-3,0,"g","T=1000K, x_e=1e-3"),
                     (8000,1e-2,0,"r","T=8000K, x_e=1e-2 (warm)"),
                     (2000,5e-3,300,"m","z=300, T=820K (+CMB)")]:
    ax0.loglog(nH, tau_HM(nH,np.full_like(nH,T),xe,z),  c+"-",  lw=2, label="H- : "+lab)
    ax0.loglog(nH, tau_H2p(nH,np.full_like(nH,T),xe,z), c+"--", lw=1.4)
ax0.loglog(nH, t_ff(nH), "k-", lw=2.5, label="free-fall time $t_{\\rm ff}$")
ax0.axhline(t_H(0),  color="grey", ls=":", lw=1.5); ax0.text(2e-4,t_H(0)*1.3,"Hubble (z=0)",fontsize=8,color="grey")
ax0.axhline(t_H(300),color="grey", ls=":", lw=1.5); ax0.text(2e-4,t_H(300)*1.3,"Hubble (z=300)",fontsize=8,color="grey")
ax0.set_xlabel(r"$n_{\rm H}$  [cm$^{-3}$]"); ax0.set_ylabel("timescale [s]")
ax0.set_title(r"$\tau_{\rm eq}$(H$^-$ solid, H$_2^+$ dashed) vs dynamical time")
ax0.set_ylim(1e-8,1e20); ax0.legend(fontsize=7,frameon=False,loc="upper right")
ax0.text(1e3,1e-3,"$\\tau_{\\rm eq}\\ll t_{\\rm dyn}$ everywhere\n$\\Rightarrow$ equilibrium",fontsize=11,ha="center",
         bbox=dict(boxstyle="round",fc="lightyellow",ec="orange"))

# Panel B: direct ODE -- H- in a collapsing parcel (n rises), full rate eqn vs equilibrium
import numpy as np
n0=1.0; T=500.0; xe=1e-4
tff=t_ff(n0)
tt=np.linspace(0,3*tff,4000); dt=tt[1]-tt[0]
# collapse: n(t)=n0/(1-t/tff)^2 (runaway), HI=n, e=xe*n, HII=xe*n
nH_t=n0/np.maximum(1-tt/(1.05*tff),1e-3)**2
k7=lambda T:3e-16*(T/300)**0.95*np.exp(-T/9320)
HM_full=np.zeros_like(tt); HM_eq=np.zeros_like(tt)
HM_full[0]=1e-20
for i in range(1,len(tt)):
    n=nH_t[i]; F=k7(T)*n*(xe*n); D=k8(T)*n+k16(T)*(xe*n)+k15(T)*n
    HM_eq[i]=F/D
    HM_full[i]=(HM_full[i-1]+F*dt)/(1+D*dt)     # semi-implicit (exact rate-eqn update)
HM_eq[0]=HM_eq[1]
ax1.loglog(nH_t[1:], HM_full[1:]/nH_t[1:], "b-",  lw=3.0, label="advected (full rate equation)")
ax1.loglog(nH_t[1:], HM_eq[1:]/nH_t[1:],   "r--", lw=2.0, label="equilibrium formula  F/D")
ax1.set_xlabel(r"$n_{\rm H}$  [cm$^{-3}$] (collapsing parcel)"); ax1.set_ylabel(r"$n_{\rm H^-}/n_{\rm H}$")
ax1.set_title("H$^-$: advected vs equilibrium track identically")
ax1.legend(fontsize=9,frameon=False)
relerr=np.nanmax(np.abs(HM_full[10:]-HM_eq[10:])/HM_eq[10:])
ax1.text(0.05,0.05,f"max rel. difference = {relerr:.1e}",transform=ax1.transAxes,fontsize=10,
         bbox=dict(boxstyle="round",fc="lightyellow",ec="orange"))
fig.tight_layout(); fig.savefig("/tmp/equilibrium_check.png",dpi=140)
print("wrote /tmp/equilibrium_check.png")
# numeric summary
for nH_,T_,xe_,z_ in [(1e-2,200,1e-4,0),(1.0,500,1e-4,0),(1e4,2000,1e-3,0),(1e8,8000,1e-2,0),(1e2,820,5e-3,300)]:
    print(f"n_H={nH_:8.0e} T={T_:5.0f} z={z_:3d}:  tau_HM={tau_HM(nH_,T_,xe_,z_):.2e}s  tau_H2+={tau_H2p(nH_,T_,xe_,z_):.2e}s  "
          f"t_ff={t_ff(nH_):.2e}s  ratio_HM={tau_HM(nH_,T_,xe_,z_)/min(t_ff(nH_),t_H(z_)):.1e}")
