# Literature notes (radiation of the runaway pair) — working file

Accessed (via WebFetch text extraction; direct PDF download blocked by proxy 403):
- Pravda & Pravdová, "Uniformly accelerated sources in electromagnetism and gravity", arXiv:gr-qc/9806114 (1998).
  * General BRS news (citing Bičák & Pravdová, gr-qc/9808068): c_{,u} = K(w)/u², w = sinθ/u.
  * Bonnor–Swaminarayan (BS) metric eq.(18)-(19): R = (ρ²+z²−t²)/2, R_i = sqrt((R−h_i)² + 2ρ²h_i);
    eq.(20) a1=(h1−h2)²/(2h2), a2=−(h1−h2)²/(2h1), k=1  ->  "two pairs of freely moving positive and
    negative particles" (no nodal singularities / struts).
  * Small-mass limit h2=h, h1=h+ε: m(1)=ε²/(2h sqrt(2h+2ε)) >0, m(2)=−ε²/((2h+2ε) sqrt(2h)) <0.
  * News: c_{,u} = ε³ 3u sin²θ / (2h (u²+2h sin²θ)^{5/2});  K(w)=±ε³ 3w²/(2h(1+2hw²)^{5/2}) (+ for u>0).
  * Bondi mass (o = 3ε³/2h), u>0:
    m = o² u³(5u⁴+32u²h+64h²) √2 (−ln u² + ln|−2√(2h)√(u²+2h)+u²+4h|) / (4096 h^{7/2}(u²+2h)^{7/2})
        + o²(u²+4h)(15u⁴+16u²h+32h²)/(3072 h³ u (u²+2h)³)
    (text extraction may be garbled; verified independently in bondi_mass_check.py)
  * "m as a function of u is everywhere decreasing. The system does radiate gravitational waves."
  * EM analogue: radiated rate constant in u* (boost symmetry); R(0,0)=2/3 (Larmor).
- Bičák & Kofroň arXiv:1210.5870 abstract only (no Bondi-mass content).
- Bičák, "Selected solutions...", gr-qc/0004016: only intro/outline of Sec. 11 extracted (BRS spacetimes are
  "the only known radiative solutions describing finite sources", admit global smooth sections of scri).
Not accessed: Bičák 1968 Proc R Soc A 302 201; Bičák–Hoenselaers–Schmidt 1983; Bičák–Schmidt 1989; Bonnor–Swaminarayan 1964.

## Progress log
- bondi_mass_check.py (done): PP closed form m(u) verified: dm/du = -(1/2)∫c_u² sinθdθ to 20 digits;
  m(u) = ∫_u^∞ L du' (so m(+∞)=0), m(-u) = -m(u), u·m → 3/8 (D=a=1) as u→0+, m ~ u^-7 at large u.
  Mass identification of the BS ± pair: A = m/d² exactly, D=md=ε³/(2h)². Turnaround-emission luminosity
  = 24/35 D²A⁴ = Phase-3 octupole value.
- linearized_news.py (partial output): linearized single hyperbolic ± dipole news = −BS news (exactly, at
  random points); the mirror pair doubles it -> ratio (both pairs)/BS = −2. The factor-2 normalisation
  mismatch is unresolved (likely the mass identification of PP's small-ε masses or c ↔ h convention).
