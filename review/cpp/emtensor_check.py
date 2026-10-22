# Independent numeric check of BiComplexScalarField::compute_emtensor
# (transcribed formulas) against the covariant stress tensor
#   T_mn = sum_s eps_s [ sum_k d_m phi_k d_n phi_k - g_mn (1/2 sum_k (dphi_k)^2 + V(|Phi_s|^2)) ]
# projected on the Eulerian observer n^mu = (1/alpha, -beta^i/alpha).
# Pi is taken as the code defines it: d_t phi = alpha Pi + beta^i d_i phi.
import numpy as np
rng = np.random.default_rng(1)
m, lam, mu = 1.0, 10240.0, 21845333.0
def V(m2): return 0.5*m*m*m2 - 0.25*lam*m2**2 + mu/6*m2**3
# random 3-metric gamma = h/chi with det h = 1
A = rng.normal(size=(3,3))*0.2; h = np.eye(3) + A@A.T*0.3
h /= np.linalg.det(h)**(1/3); chi = 0.8
gam = h/chi; gamU = np.linalg.inv(gam); hU = np.linalg.inv(h)
alpha = 0.93; beta = rng.normal(size=3)*0.1
g = np.zeros((4,4)); betaL = gam@beta
g[0,0] = -alpha**2 + beta@betaL; g[0,1:] = g[1:,0] = betaL; g[1:,1:] = gam
gU = np.linalg.inv(g)
n_up = np.array([1/alpha, *(-beta/alpha)])
sectors = []
for s in (+1, -1):
    phi = rng.normal(size=2)*0.02; Pi = rng.normal(size=2)*0.02
    dphi = rng.normal(size=(2,3))*0.01
    sectors.append((s, phi, Pi, dphi))
# covariant
T = np.zeros((4,4))
for s, phi, Pi, dphi in sectors:
    for k in range(2):
        dt = alpha*Pi[k] + beta@dphi[k]
        d4 = np.array([dt, *dphi[k]])
        T += s*np.outer(d4, d4)
        T -= s*g*0.5*(d4@gU@d4)
    T -= s*g*V(phi@phi)
rho_c = n_up@T@n_up
S_c = -(T[1:,:]@n_up)        # S_i = -gamma_i^mu n^nu T_mu nu (spatial lower i)
Sij_c = T[1:,1:]
trS_c = np.einsum('ij,ij', gamU, Sij_c)
# code transcription
rho = 0; j = np.zeros(3); S = np.zeros((3,3))
for s, phi, Pi, dphi in sectors:
    for k in range(2):
        Vt = -Pi[k]**2 + chi*dphi[k]@hU@dphi[k]
        rho += s*(Pi[k]**2 + 0.5*Vt)
        j += s*(-Pi[k]*dphi[k])
        S += s*(-0.5*h*Vt/chi + np.outer(dphi[k], dphi[k]))
    rho += s*V(phi@phi); S -= s*h*V(phi@phi)/chi
trS = chi*np.einsum('ij,ij', hU, S)
print("rho  ", rho, rho_c, abs(rho-rho_c))
print("S_i  ", np.max(abs(j-S_c)))
print("S_ij ", np.max(abs(S-Sij_c)))
print("trS  ", trS, trS_c, abs(trS-trS_c))
# potential derivative check: dV/dphi1 vs code 2 phi1 dV/d|Phi|^2
p = np.array([0.013, -0.007]); eps=1e-7
num = (V((p+[eps,0])@(p+[eps,0])) - V((p-[eps,0])@(p-[eps,0])))/(2*eps)
m2 = p@p; code = 2*p[0]*(0.5*m*m - 0.5*lam*m2 + 0.5*mu*m2*m2)
ode = m*m*p[0] - lam*m2*p[0] + mu*m2*m2*p[0]
print("dV/dphi1 numeric", num, "code", code, "ID-ODE form", ode)
