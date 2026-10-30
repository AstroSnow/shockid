"""3D validation of shockid against exact MHD Rankine-Hugoniot shocks.

Each shock is solved exactly (coplanar, gamma=5/3, B in code units with B^2/2 pressure),
smoothed with a tanh profile, rotated to an arbitrary normal in 3D, given a lab-frame
shock speed, and run through shockid. The dominant detected class and the median
measured Mach numbers are compared with the exact values.
"""
import numpy as np
from scipy.optimize import brentq
from shockid import shockid, getState, TRANS

g = 5.0/3.0

def fluxes(s, Bn):
    ro, u, vt, bt, p = s
    return np.array([ro*u, ro*u*u + p + bt*bt/2, ro*u*vt - Bn*bt, u*bt - Bn*vt,
                     u*(g/(g-1)*p + 0.5*ro*(u*u+vt*vt) + bt*bt + Bn*Bn) - Bn*(u*Bn + vt*bt)])

def downstream_states(up, Bn):
    """All compressive RH solutions for upstream (ro,u,vt,bt,p) in the shock frame."""
    F1 = fluxes(up, Bn); m = F1[0]
    def resid(r):
        ro = up[0]*r; u = m/ro
        vt, bt = np.linalg.solve([[m, -Bn], [-Bn, u]], [F1[2], F1[3]])
        p = F1[1] - m*u - bt*bt/2
        s = np.array([ro, u, vt, bt, p])
        return fluxes(s, Bn)[4] - F1[4], s
    rs = np.linspace(1.001, 3.99, 4000); e = np.array([resid(r)[0] for r in rs])
    out = []
    for i in np.where(np.sign(e[:-1]) != np.sign(e[1:]))[0]:
        s = resid(brentq(lambda r: resid(r)[0], rs[i], rs[i+1]))[1]
        if s[4] > 0: out.append(s)
    return out

def mach(s, Bn):
    ro, u, vt, bt, p = s; cs2 = g*p/ro; b2 = Bn**2 + bt**2; va2 = b2/ro
    d = np.sqrt((cs2+va2)**2 - 4*va2*cs2*Bn**2/b2)
    return np.array([abs(u)/np.sqrt(0.5*(cs2+va2-d)), abs(u)/np.sqrt(Bn**2/ro),
                     abs(u)/np.sqrt(0.5*(cs2+va2+d))])

def embed(up, dn, Bn, n, t, nx, ny, nz, Vs=0.3, w=0.03):
    xg, yg, zg = (np.linspace(0, 1, k) for k in (nx, ny, nz))
    Z, Y, X = np.meshgrid(zg, yg, xg, indexing='ij')
    f = 0.5*(1 + np.tanh(((X-.5)*n[0] + (Y-.5)*n[1] + (Z-.5)*n[2])/w))
    ro, u, vt, bt, p = [a + (b-a)*f for a, b in zip(up, dn)]
    u = u + Vs
    v = [u*n[i] + vt*t[i] for i in range(3)]
    B = [Bn*n[i] + bt*t[i] for i in range(3)]
    return xg, yg, zg, [ro, *v, *B, p]

if __name__ == '__main__':
    n = np.array([1, 0.6, 0.8]); n /= np.linalg.norm(n)
    t = np.cross(n, [0, 0, 1.]); t /= np.linalg.norm(t)
    cases = [('fast',         np.array([1, 3.0, 0, 0.8, 0.6]), 0.6),
             ('slow',         np.array([1, 0.9, 0, 0.6, 0.15]), 1.0),
             ('slow',         np.array([1, 0.8, 0, 0.5, 0.3]), 1.2),
             # 1->3 needs low beta to sit clearly away from M_A = 1 (post M_A = 0.85)
             ('intermediate', np.array([1, 1.39, 0, 0.25, 0.02]), 1.0),
             ('intermediate', np.array([1, 1.2, 0, 0.05, 0.05]), 1.0)]
    options = [{}, {'dbproject': True}]
    for label, up, Bn in cases:
        for dn in downstream_states(up, Bn):
            if label == 'intermediate' and dn[3]*up[3] >= 0:
                continue                      # keep only the Bt-reversing branch
            if np.min(np.abs(np.r_[mach(up, Bn), mach(dn, Bn)] - 1)) < 0.05:
                continue                      # marginal by construction: getState rejects it
            mp, mq = mach(up, Bn), mach(dn, Bn)
            want = TRANS.get((getState(*mp), getState(*mq)))
            for grid, kw in [(g_, o_) for g_ in [(24, 24, 24), (16, 32, 24)] for o_ in options]:
                xg, yg, zg, d = embed(up, dn, Bn, n, t, *grid)
                sh = shockid(xg, yg, zg, *d, ndim=3, convl=1e-3, avecyl=2, **kw)
                cnt = {k: len(sh[k]) for k in TRANS.values()}
                got = max(cnt, key=cnt.get)
                frac = cnt[got]/max(sum(cnt.values()), 1)
                print(f'{label:12s} grid={grid} {str(kw):20s} expected={want:5s} got={got:5s} '
                      f'({100*frac:.0f}% of hits)  exact fast Mach {mp[2]:.2f}->{mq[2]:.2f}, '
                      f'measured {np.median(sh[got+"mach"][:,2]):.2f}->{np.median(sh[got+"mach"][:,5]):.2f}'
                      + ('' if got == want else '   <-- MISMATCH'))
