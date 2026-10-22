"""Effect of MU = 21845333 (paper literal) vs LAM^2/4.8 = 21845333.33 (ours) on the canonical star."""
import iso_star as I
for mu in (I.MU_PAPER, I.MU_EXACT):
    e = I.at_omega(0.750, +1, mu=mu, which='w_true')
    print(f"mu={mu:.4f}: phic={e['phic']:.10f} M_true={e['M_true']:+.8f} M_paper-style={e['M_paper']:+.8f}")
