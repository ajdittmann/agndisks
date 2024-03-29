import numpy as np
import argparse as ag
import matplotlib
matplotlib.rc('text', usetex=True)
import matplotlib.pyplot as plt
from numba import njit  #improves things if Nr > 10^5 or so
from scipy import optimize as sciop
from scipy import interpolate as scint

# resolution
rmax = 10**6.0	#maximum radius in Rg
rmin = 6	#minimum radius in Rg
Nr = 10000

#model parameters
eta = 0.1	# Eddington fraction
logm = 7	# log_10(MBH/Msun)
alpha = 0.1	# effective viscosity parameter
X = 0.72	# hydrogen mass fraction
Z = 0.02	# 'metals' mass fraction
eps = 0.1	# efficiency of rest mass -> energy for eddington accretion rate

#opacity table - set to None for Kramers opacity
opacTable = None
#opacTable = "combined-opacity/opacitysolar09dustq3p5amax0p1new.txt"


#physical constants
MSUN = 1.989*10**33
G = 6.6743*10**-8
c = 2.99792458*10**10
kb = .3807*10**-16
sigma = 5.6704*10**-5
mh = 1.6726*10**-24
sigmaT = 6.6524587158*10**-25
a = 4*sigma/c

#Derived parameters
M = 10.0**logm
MBH = M*MSUN
Mdot = eta*4*np.pi*G*MBH*mh/(eps*c*sigmaT)
mu = mh*4/(3+5*X-Z)
Rg = G*MBH/c**2
rin = 3*Rg
kes = 0.2*(1 + X)

@njit
def Omega(r, M):
  return np.sqrt(G*M/r**3)

@njit
def opacKramer(rho, T):
  kk = 4.0*10**25*(X+1)*(Z + 0.001)*rho*T**-3.5
  return kes + kk

if opacTable is not None:
  with open(opacTable,'r') as fp:
    line=fp.readline()
    numbers_str = line.split()
    nrho=int(numbers_str[1])
    nt=int(numbers_str[3])
    rhoread=np.zeros(nrho)
    treadnew=np.zeros(nt)
    rosscombineread=np.zeros([nt,nrho])
    planckcombineread=np.zeros([nt,nrho])
    for irho in range(nrho):
        for it in range(nt):
            line=fp.readline()
            numbers_str = line.split()
            rhoread[irho]=float(numbers_str[0])
            treadnew[it]=float(numbers_str[1])
            rosscombineread[it,irho]=float(numbers_str[2])
  opacityinter = scint.RectBivariateSpline(np.log10(rhoread), np.log10(treadnew), np.log10(rosscombineread.T), kx=1, ky=1)
  maxT = np.max(treadnew)

  def opacTabulated(rho, T):
    if T > maxT: return kes
    else: return 10.0**opacityinter.ev(np.log10(rho), np.log10(T))

def opacity(rho, T):
  if opacTable is not None:
    opac = opacTabulated(rho, T)
  else:
    opac = opacKramer(rho, T)
  return opac

def rhosolve(rho_g, tau, T, B):
  kappa = opacity(rho_g, T)
  return tau - kappa*B**(1./3.)*rho_g**(2./3.)

@njit
def psolve(T, rho, P):
  p = rho*kb*T/mu + a*T**4/3
  return p-P

def tsolve(Tg, A, B, omega, rho_g):
  tau = Tg**4/A
  res = sciop.fsolve(rhosolve, x0=rho_g, args = (tau, Tg, B), full_output=1)
  rho = res[0]
  H = (B/rho)**(1./3.)
  P = rho*H**2*omega**2
  res = sciop.fsolve(psolve, x0 = Tg, args=(rho,P), full_output=1)
  return res[0] - Tg

def getvals(r, T_g = None, rho_g=None):
  omega = Omega(r, MBH)
  mdot = Mdot*(1-np.sqrt(rin/r))
  A = 9*mdot*omega**2/(np.pi*64*sigma)
  B = mdot/(6*np.pi*alpha*omega)
  if T_g is None:
    T_g = omega*(3/a)**0.5*(B/(A*kes))**0.25
  if rho_g is None:
    rho_g = omega**6*B*(3/(a*A*kes))**3
  resT = sciop.fsolve(tsolve, x0=T_g, args=(A, B, omega, rho_g), full_output=1)
  T = resT[0]
  tierr = resT[2]
  tau = T**4/A
  resrho = sciop.fsolve(rhosolve, x0=rho_g, args = (tau, T_g, B), full_output=1 )
  rho = resrho[0]
  rhoierr = resrho[2]
  H = (B/rho)**(1./3.)
  P = rho*H**2*omega**2
  return T, rho, H, P, tau, tierr, rhoierr

rs = 10.0**np.linspace(np.log10(rmin), np.log10(rmax), Nr)*Rg
ts = np.empty(Nr)
rhos = np.empty(Nr)
hs = np.empty(Nr)
ps = np.empty(Nr)
taus = np.empty(Nr)
badTs = np.empty(Nr)
badRhos = np.empty(Nr)
for i in range(Nr):
  if i==0:
    Ti, rhoi, Hi, Pi, taui, terr, rerr = getvals(rs[i])
  else:
    if badTs[i-1]==1: tguess = ts[i-1]
    else:
      lastGoodT = np.nonzero(badTs==1)[0][-1]
      tguess = ts[lastGoodT]
    if badRhos[i-1]==1: rguess = rhos[i-1]
    else:
      lastGoodR = np.nonzero(badRhos==1)[0][-1]
      rguess = rhos[lastGoodR]
    Ti, rhoi, Hi, Pi, taui, terr, rerr = getvals(rs[i], tguess, rguess)

  hs[i] = Hi
  taus[i] = taui
  ps[i] = Pi
  badTs[i] = terr
  badRhos[i] = rerr
  ts[i] = Ti
  rhos[i] = rhoi

#only use radii where it actually found a valid solution
good = (badTs==1)&(badRhos==1)

rs = rs[good]
ts = ts[good]
rhos = rhos[good]
hs = hs[good]
taus = taus[good]
ps = ps[good]
Nr = np.sum(good)

omegas = Omega(rs, MBH)
cs = hs*omegas
kappas = taus/(rhos*hs)
sigmas = 2*hs*rhos
Qs = cs*omegas/(np.pi*sigmas*G)

# below are some example plotting functions for the disk model outputs.
# for setting Rg and pc axis labels simultaneously
def fwd(rg):
  return rg*Rg/(3.086*10**18)
def bwd(pc):
  return pc*(3.086*10**18)/Rg

fig, ax = plt.subplots(4,2, figsize=(5.0, 5.0))
ax[0,0].loglog(rs/Rg, rhos)
ax[1,0].loglog(rs/Rg, ts)
ax[2,0].loglog(rs/Rg, sigmas)
ax[3,0].loglog(rs/Rg, cs)


ax00pc = ax[0,0].secondary_xaxis('top', functions=(fwd, bwd))

ax[0,1].loglog(rs/Rg, kappas)
ax[1,1].loglog(rs/Rg, taus)
ax[1,1].loglog(rs/Rg, np.ones(Nr), color='black', ls = "--")
ax[2,1].loglog(rs/Rg, hs/rs)
ax[3,1].loglog(rs/Rg, Qs)
ax[3,1].loglog(rs/Rg, np.ones(Nr), color='black', ls = "--")

ax01pc = ax[0,1].secondary_xaxis('top', functions=(fwd, bwd))


plt.subplots_adjust(top=0.9125, bottom=0.065, hspace=0.0, wspace=0.3, right =0.99)

ax[3,0].set_xlabel(r'$r/r_g$',labelpad=-3)
ax[3,1].set_xlabel(r'$r/r_g$',labelpad=-3)
ax01pc.set_xlabel(r'$r~[pc]$')
ax00pc.set_xlabel(r'$r~[pc]$')

ax[0,0].set_ylabel(r'$\rho$')
ax[1,0].set_ylabel(r'$T$')
ax[2,0].set_ylabel(r'$\Sigma$')
ax[3,0].set_ylabel(r'$c_s$')

ax[0,1].set_ylabel(r'$\kappa$')
ax[1,1].set_ylabel(r'$\tau$')
ax[2,1].set_ylabel(r'$H/r$',labelpad=-2)
ax[3,1].set_ylabel(r'$Q$')

plt.savefig('diskmodel.pdf')

output = np.empty((Nr,12))
output[:,0] = rs/Rg
output[:,1] = fwd(rs/Rg)
output[:,2] = rhos
output[:,3] = omegas
output[:,4] = ts
output[:,5] = ps
output[:,6] = sigmas
output[:,7] = cs
output[:,8] = kappas
output[:,9] = taus
output[:,10] = hs/rs
output[:,11] = Qs
head = "0:r/Rg\t1:r/pc\t2:density\t3:angular velocity\t4:temperature\t5:pressure\t6:surface density\t7:sound speed\t8:opacity\t9:optical depth\t10:disk aspect ratio (H/r)\t11:Toomre Q"
np.savetxt('diskmodel.txt', output, header=head)
