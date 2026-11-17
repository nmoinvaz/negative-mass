# winding2 normalisation on a box that contains a KNOWN zero: the mini-BS QNM reflected... QNMs are in Re s<0,
# so instead test the canonical-star crossing: radial agent's unstable canonical l=0 is not polar l>=2.
# Use the mini star with the box shifted to contain s_qnm (Re s<0 allowed here since kp,km cuts are on Re s=0
# only for |Im s| > 1-w ... box is small and away from them).
#
# Options (the defaults run the two boxes of the first version):
#   --strip         instead, the strip geometry of the 0 < Re s < 0.1 counts: boxes -2e-4 < Re s < 0.1 below the
#                   threshold Im s = 1-w = 0.1469 that hold the upper QNM conj(s_qnm) and the zeros above it that
#                   polish_zero.py --mini found (see the list in the code): 1, 3, 2 and 0 zeros
#   --rsmult=40     R_s = r_m + rsmult / (min Re k on the box)
#   --cachedir=DIR  keep the evaluations of D in DIR (restart-safe; default: none)
#   --mutate=NAME   reverse (traverse clockwise) or shift (move the first box off the zero, expectation unchanged);
#                   the script must exit 1
# Exit status 1 unless every winding equals the number of known zeros in its box, to within 0.05.
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pickle, mini_bs as mb, background as B
from polar_solver import PolarProblem
opts = dict((a[2:].split('=', 1) + ['1'])[:2] for a in sys.argv[1:] if a.startswith('--'))
assert set(opts) <= {'strip', 'rsmult', 'mutate', 'cachedir'}, opts
strip, rsmult, mutate = 'strip' in opts, float(opts.get('rsmult', 40)), opts.get('mutate')
assert mutate in (None, 'reverse', 'shift'), mutate
st=mb.star_from_shoot(0.08)
for _ in range(3): st=B.solve_star(0.853087,+1,guess=st,R=st.R)
P=PolarProblem(st,2,LAM=0.,MU=0.); w=st.w
c=-5.353e-5-0.1195191j; hx,hy=2e-4,1.5e-3
ks=lambda s:(np.sqrt(1+(s-1j*w)**2),np.sqrt(1+(s+1j*w)**2))
def Dc_raw(s):
    kp,km=ks(s); return P.D(complex(s))*np.exp(-1j*((kp+km)*(P.Rs_fixed-P.rm)+(kp+km+s)*P.rm).imag)
cache={}
def Dc(s):
    s=complex(round(s.real,12),round(s.imag,12))
    if s not in cache:
        cache[s]=Dc_raw(s)
        if cf: pickle.dump(cache,open(cf,'wb'))
    return cache[s]
cf=None
if not strip:
    boxes=[(c,'around QNM',1),(c+0.004j,'displaced',0)]
    corners=lambda cc:[cc-hx-1j*hy,cc+hx-1j*hy,cc+hx+1j*hy,cc-hx+1j*hy,cc-hx-1j*hy]
    cases=[(lab,corners(cc+(0.02j if mutate=='shift' and lab=='around QNM' else 0)),ex) for cc,lab,ex in boxes]
else:
    # zeros of D in Im s in (0.09, 0.15) near the axis, located by polish_zero.py --mini (secant, rtol 1e-8 and 1e-10):
    #   -5.35350e-5 + 0.11951908 i  (the published QNM), -2.62464e-5 + 0.131605 i, -1.35809e-5 + 0.137308 i,
    #   -7.7455e-6 + 0.140372 i  (a series that accumulates at the threshold Im s = 1 - w = 0.1469)
    xa,xb=-2e-4,0.1
    def rect(ya,yb): return [xa+1j*ya,xb+1j*ya,xb+1j*yb,xa+1j*yb,xa+1j*ya]
    dy=-0.05 if mutate=='shift' else 0
    cases=[('strip, fundamental only',rect(0.09+dy,0.125+dy),1),
           ('strip, three zeros',rect(0.09+dy,0.14+dy),3),
           ('strip, above the QNM',rect(0.1225,0.14),2),
           ('strip, between the QNM and the next zero',rect(0.1225,0.1285),0)]
if mutate=='reverse': cases=[(lab,cor[::-1],ex) for lab,cor,ex in cases]
def box_winding(cor, nseg=8):
    """winding of Dc around the polygon cor; a step is accepted only if its two half steps are small and add up to it
    (for the first two bisection levels), see winding2.py"""
    n=[0]; mx=[0.0]
    def step(a,b,t0,t1,depth):
        da,db=Dc(a+(b-a)*t0),Dc(a+(b-a)*t1); d=np.angle(db/da); n[0]+=1
        fine=t1-t0<=1e-5; tm=0.5*(t0+t1)
        if abs(d)<=0.4 and depth<=1 and not fine:
            dm=Dc(a+(b-a)*tm); h1,h2=np.angle(dm/da),np.angle(db/dm)
            if max(abs(h1),abs(h2))>0.4 or abs(h1+h2-d)>0.1: return step(a,b,t0,tm,depth+1)+step(a,b,tm,t1,depth+1)
        if abs(d)>0.4 and not fine: return step(a,b,t0,tm,depth+1)+step(a,b,tm,t1,depth+1)
        mx[0]=max(mx[0],abs(d)); return d
    tot=0.0
    for a,b in zip(cor[:-1],cor[1:]):
        tt=np.linspace(0,1,nseg+1)
        for t0,t1 in zip(tt[:-1],tt[1:]): tot+=step(a,b,t0,t1,0)
    return tot/2/np.pi,n[0],mx[0]
ok=True
for lab,cor,ex in cases:
    kmin=min(min(k.real for k in ks(cor[i]+(cor[i+1]-cor[i])*t)) for i in range(4) for t in np.linspace(0,1,200)) if strip \
        else min(k.real for k in ks(c))
    P.Rs_fixed=P.rm+rsmult/kmin
    cache={}; cf=None
    if 'cachedir' in opts:
        cf=os.path.join(opts['cachedir'],f"control_cache_{'strip' if strip else 'box'}_{lab.replace(' ','_').replace(',','')}_Rs{P.Rs_fixed:.1f}.pkl")
        if os.path.exists(cf): cache=pickle.load(open(cf,'rb'))
    wd,n,mx=box_winding(cor); good=abs(wd-ex)<0.05; ok&=good
    print(f"winding2 normalisation control, {lab}: winding = {wd:+.4f}  (expected {ex}, R_s = {P.Rs_fixed:.0f}, "
          f"min Re k {kmin:.3f}, rsmult {rsmult:g}, {n} steps, max step {mx:.2f} rad)  {'PASS' if good else 'FAIL'}"
          +(f"  [mutate={mutate}]" if mutate else ""),flush=True)
sys.exit(0 if ok else 1)
