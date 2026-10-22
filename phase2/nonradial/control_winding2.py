# winding2 normalisation on a box that contains a KNOWN zero: the mini-BS QNM reflected... QNMs are in Re s<0,
# so instead test the canonical-star crossing: radial agent's unstable canonical l=0 is not polar l>=2.
# Use the mini star with the box shifted to contain s_qnm (Re s<0 allowed here since kp,km cuts are on Re s=0
# only for |Im s| > 1-w ... box is small and away from them).
import sys; sys.path.insert(0,'/home/claude/negative-mass/phase2/nonradial')
import numpy as np, mini_bs as mb, background as B
from polar_solver import PolarProblem
st=mb.star_from_shoot(0.08)
for _ in range(3): st=B.solve_star(0.853087,+1,guess=st,R=st.R)
P=PolarProblem(st,2,LAM=0.,MU=0.); w=st.w
c=-5.353e-5-0.1195191j; hx,hy=2e-4,1.5e-3
ks=lambda s:(np.sqrt(1+(s-1j*w)**2),np.sqrt(1+(s+1j*w)**2))
P.Rs_fixed=P.rm+40/min(k.real for k in ks(c))
def Dc(s):
    kp,km=ks(s); return P.D(complex(s))*np.exp(-1j*((kp+km)*(P.Rs_fixed-P.rm)+(kp+km+s)*P.rm).imag)
for cc,lab in [(c,'around QNM'),(c+0.004j,'displaced')]:
    cor=[cc-hx-1j*hy,cc+hx-1j*hy,cc+hx+1j*hy,cc-hx+1j*hy,cc-hx-1j*hy]; tot=0; n=0
    for a,b in zip(cor[:-1],cor[1:]):
        seg=list(np.linspace(0,1,9)); i=0
        while i<len(seg)-1:
            d=np.angle(Dc(a+(b-a)*seg[i+1])/Dc(a+(b-a)*seg[i])); n+=1
            if abs(d)>0.4: seg.insert(i+1,0.5*(seg[i]+seg[i+1])); continue
            tot+=d; i+=1
    print(f"winding2 normalisation control, {lab}: winding = {tot/2/np.pi:+.4f}")
