"""
Order-of-magnitude radiation estimate for the uniformly accelerated Bondi dipole (Newtonian multipoles,
quadrupole formalism of Thorne 1980, G = c = 1). Bodies on the z axis: z1 = Z + d/2 (mu = +m), z2 = Z - d/2
(mu = -m), Z = a t^2/2. Mass multipoles I_L = sum mu x^<L> (STF). Current moments vanish (x x v = 0).
Luminosities: L2 = (1/5) I3dot_ij I3dot_ij,  L3 = (1/189) I4dot_ijk I4dot_ijk.
"""
import sympy as sp
import itertools as it
t, a, d, m = sp.symbols('t a d m', positive=True)
Z = a * t**2 / 2
bodies = ((m, Z + d / 2), (-m, Z - d / 2))
def stf2(z):
    x = [0, 0, z]; r2 = z * z
    return [[x[i] * x[j] - sp.Rational(1, 3) * int(i == j) * r2 for j in range(3)] for i in range(3)]
def stf3(z):
    x = [0, 0, z]; r2 = z * z
    return {(i, j, k): x[i] * x[j] * x[k] - sp.Rational(1, 5) * r2 * (int(i == j) * x[k] + int(i == k) * x[j] + int(j == k) * x[i])
            for i, j, k in it.product(range(3), repeat=3)}
I2 = [[sum(mu * stf2(z)[i][j] for mu, z in bodies) for j in range(3)] for i in range(3)]
I3 = {key: sum(mu * stf3(z)[key] for mu, z in bodies) for key in it.product(range(3), repeat=3)}
L2 = sp.simplify(sum(sp.diff(I2[i][j], t, 3)**2 for i in range(3) for j in range(3)) / 5)
L3 = sp.simplify(sum(sp.diff(v, t, 4)**2 for v in I3.values()) / 189)
print("dipole moment 2nd derivative:", sp.simplify(sp.diff(sum(mu * z for mu, z in bodies), t, 2)))
print("quadrupole luminosity L2 =", L2)
print("octupole  luminosity L3 =", L3, " = (with a = m/d^2)", sp.simplify(L3.subs(a, m / d**2)))
num = {m: 0.0143, d: 10}
L3n = float(L3.subs(a, m / d**2).subs(num))
print(f"m = 0.0143, d = 10: L3 = {L3n:.2e};  energy radiated by t = 1000: {L3n*1000:.1e}  vs  E_tot = m^2/d = {0.0143**2/10:.1e}")
print("(The slow-motion expansion about a fixed origin fails once v ~ 1; higher multipoles grow like powers of v.)")
