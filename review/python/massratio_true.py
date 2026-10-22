"""True vs 5%-radius masses for the mass-ratio cells' retuned stars (trusted iso solver)."""
import sys
sys.path.insert(0, '/home/claude/negative-mass/phase2/mass_check')
import iso_star as I
csv = {('p',0.804):-0.011472499, ('c',0.81):0.010720519, ('p',0.88):-0.008573106}
MT_P, MT_M, M_P, M_M = 0.0146265, 0.0145803, 0.014349964, 0.014294692
res = {}
for (s, w), mc in csv.items():
    e = I.at_omega(w, +1 if s=='c' else -1, which='w_paper'); res[(s,w)] = e['M_true']
    print(f"{s} w={w}: csv {mc:+.7f} ours-5% {e['M_paper']:+.7f} true {e['M_true']:+.7f} csv/true {mc/e['M_true']:.5f}", flush=True)
pred_paper = {'w0804': 0.011472499/M_M, 'w088': 0.008573106/M_M, 'heavy': 0.010720519/M_P}
pred_true = {'w0804': -res[('p',0.804)]/MT_M, 'w088': -res[('p',0.88)]/MT_M, 'heavy': res[('c',0.81)]/MT_P}
for k in pred_paper: print(f"{k}: predicted pull ratio paper-M {pred_paper[k]:.4f}  true-M {pred_true[k]:.4f}  shift {pred_true[k]/pred_paper[k]-1:+.2%}")
rat_p = {'heavy': M_M/0.010720519, 'w0804': 0.011472499/M_P, 'w088': 0.008573106/M_P}
rat_t = {'heavy': MT_M/res[('c',0.81)], 'w0804': -res[('p',0.804)]/MT_P, 'w088': -res[('p',0.88)]/MT_P}
for k in rat_p: print(f"{k}: |M-|/M+ paper {rat_p[k]:.4f} true {rat_t[k]:.4f}")
print(f"matched: paper {M_M/M_P:.4f} true {MT_M/MT_P:.4f}")
