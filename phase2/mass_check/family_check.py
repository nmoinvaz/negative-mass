"""Spot-check two more points of the paper's stars/star_family.csv with both extractions."""
import iso_star as I
for w_t, eps, Mp in [(0.80, +1, 0.011209378), (0.70, -1, -0.020520688)]:
    e = I.at_omega(w_t, eps, which='w_paper')
    print(f"eps={eps:+d} w={w_t}: phic={e['phic']:.9f}  M_paper-style={e['M_paper']:+.9f} (paper {Mp:+.9f}, "
          f"rel {e['M_paper']/Mp-1:+.1e})  M_true={e['M_true']:+.9f}  paper/true={Mp/e['M_true']:.5f}", flush=True)
