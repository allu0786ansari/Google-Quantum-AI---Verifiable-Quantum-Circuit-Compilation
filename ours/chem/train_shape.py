"""Train our own DF-THC factorization at a cost-optimal shape with the starter kit's trainer
(objective of the paper, eq. 48: minimise lambda_eff subject to the 20 mHa limit).
    python train_shape.py R B C epochs lr out.npz [seed.npz]
"""
import sys, time
import numpy as np
KIT = '/home/claude/data/starter/starter_kit'
sys.path.insert(0, KIT)
from dfthc.train_df_thc_sga import Initialize, Train
from starter.task1_factorization.train_df_thc_qaggle import (load_seed_params, frozen_symmetry_params,
                                                             save_solution, PACKAGED_FCIDUMP, PACKAGED_BASELINE)
R, B, C, epochs, lr, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), sys.argv[6]
seed = sys.argv[7] if len(sys.argv) > 7 else PACKAGED_BASELINE
eps = float(sys.argv[8]) if len(sys.argv) > 8 else 0.010
ref = Initialize(PACKAGED_FCIDUMP, -22140.27335)
params, _ = load_seed_params(seed, rank=R, bases=B, copies=C, orbitals=ref['orbitals'])
frozen = frozen_symmetry_params(ref)
params = {**params, **frozen}
t = time.time()
sol = Train(generating_set='Qaggle', method='DF-THC-SOS', epochs=epochs, learning_rate=lr, rho=1e-4,
            params=params, ref=ref, print_interval=max(1, epochs // 20), epsilon_dfthc=eps,
            frozen_params=frozen, lambda_floor=49.0, sos_gap_floor=3.914898, learning_rate_final=lr / 10)
v = sol['values']
print('done', round(time.time() - t), 's', {k: float(v[k]) for k in ('epsilon_fro', 'Lambda', 'sos_gap', 'lambda_eff')}, flush=True)
save_solution(out, sol['params_best'], v, ref=ref, metadata={'shape': [R, B, C], 'epochs': epochs, 'seed': seed})
