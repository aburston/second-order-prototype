import numpy as np, sys
from scipy.integrate import solve_ivp
wn = 1.0
def xeq(zm, v0): return 2*zm*v0/wn   # u=0: wn^2 x1 = 2 zm wn v0

def P(zp, zm, v0, r):
    """one return to the section {x2 = 0, x1 > xeq}: (next radius, period)"""
    xe = xeq(zm, v0)
    def f(t, y):
        w = y[1] - v0
        return [y[1], -wn**2*y[0] - 2*(zp if w > 0 else zm)*wn*w]
    def ev(t, y): return y[1]
    ev.direction = -1
    s = solve_ivp(f, (0, 30), [xe + r, 0.0], events=ev, rtol=1e-12, atol=1e-14)
    i = [k for k, t in enumerate(s.t_events[0]) if t > 1e-6]
    if not i: return np.nan, np.nan
    return s.y_events[0][i[0]][0] - xe, s.t_events[0][i[0]]

def orbit(zp, zm, v0, r, n=400, tol=1e-11):
    T = np.nan
    for _ in range(n):
        rn, T = P(zp, zm, v0, r)
        if not np.isfinite(rn) or rn > 1e9 or rn < 1e-9: return rn, T
        if abs(rn - r) < tol*max(1.0, abs(r)): return rn, T
        r = rn
    return r, T

zp, zm, v0 = 0.3, -0.1, 1.0
rs, _ = orbit(zp, zm, v0, 2.0)
print(f"fixed point r* = {rs:.9f},  period T = {P(zp,zm,v0,rs)[1]:.6f}")
print("Poincare map around it:")
for r in [rs*0.3, rs*0.9, rs, rs*1.1, rs*3.0]:
    pr, _ = P(zp, zm, v0, r)
    print(f"  r={r:9.6f} -> P(r)={pr:9.6f}   P(r)-r = {pr-r:+.3e}")
h = 1e-6
m = (P(zp,zm,v0,rs+h)[0] - P(zp,zm,v0,rs-h)[0])/(2*h)
print(f"  multiplier dP/dr = {m:.6f}  -> |m|<1: {abs(m)<1}   hyperbolic and attracting")
sys.stdout.flush()

print("\namplitude scales exactly linearly with the offset v0:")
for v in [0.25, 1.0, 4.0, 16.0]:
    r, _ = orbit(zp, zm, v, 3*v)
    print(f"  v0={v:>6}:  r* = {r:13.8f}   r*/v0 = {r/v:.9f}")
sys.stdout.flush()

print("\npersistence under parameter perturbation (structural stability):")
for d in [-0.03, -0.01, 0.0, 0.01, 0.03]:
    r, _ = orbit(zp+d, zm, v0, 2.0)
    print(f"  zeta+ = {zp+d:.3f}:  r* = {r:.6f}")
sys.stdout.flush()

print("\nexistence condition  zeta- < 0 < mean damping:")
for a, b in [(0.30,-0.10), (0.05,-0.10), (0.30, 0.10), (0.10,-0.30)]:
    r, _ = orbit(a, b, v0, 1.0)
    got = ("grows unbounded" if (not np.isfinite(r) or r > 1e6)
           else "decays to the equilibrium" if r < 1e-6 else f"limit cycle r* = {r:.4f}")
    print(f"  zeta+={a:+.2f} zeta-={b:+.2f} mean={(a+b)/2:+.3f}  ->  {got}")
sys.stdout.flush()

print("\ncontrast: v0 = 0, boundary back through the equilibrium, same zetas:")
r, seq = 1.0, []
for _ in range(4):
    r, _ = P(zp, zm, 0.0, r); seq.append(r)
print(f"  radii:  {[f'{x:.6f}' for x in seq]}")
print(f"  ratios: {[f'{seq[i+1]/seq[i]:.6f}' for i in range(3)]}   constant -> pure scaling, no isolated cycle")
