import json, sys, numpy as np
st = np.load(sys.argv[1]); shots = int(sys.argv[2])
ops = np.memmap('ops893.bin', dtype=np.int32, mode='r').reshape(-1, 5)
idx, act, a1, b1, fire = st.T.astype(np.int64)
code = ops[idx, 0]
tof = np.isin(code, [2, 5, 15])
idx, act, a1, b1, fire, code = idx[tof], act[tof], a1[tof], b1[tof], fire[tof], code[tof]
mean_T = act.sum() / shots
print(f'mean Toffoli-like per shot (sim): {mean_T:.1f}  ops {len(idx)}')
cats = {}
dead = fire == 0
cats['never fires (delete)'] = dead
ccx_like = np.isin(code, [2, 15])
a_one = ccx_like & ~dead & ((a1 == act) | (b1 == act))
cats['a control always 1 (-> CX/CZ, free)'] = a_one
for k, m in cats.items():
    print(f'{k:40s} ops {m.sum():7d}  mean Toffolis saved/shot {act[m].sum() / shots:9.1f}')
np.save('cand_' + sys.argv[1], dict(dead=idx[dead], a_one=idx[a_one]), allow_pickle=True)
