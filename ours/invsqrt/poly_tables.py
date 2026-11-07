"""Freeze the piecewise-cubic tables used by the circuit (K=7, T=29, z bits 15/14/11) and the model Y."""
import numpy as np
from poly_model import fit
r = fit(7, 3, {1: 29, 2: 29, 3: 29}, {1: 15, 2: 14, 3: 11}, 29, tables=True)
print('worst rel', np.log2(r['rel']), 'widths', r['widths'], 'Y bits', int(r['Y'].max()).bit_length(), 'Ymin', r['Y'].min())
np.savez('poly_tables.npz', C0=r['C'][0], C1=r['C'][1], C2=r['C'][2], C3=r['C'][3], Y=r['Y'],
         K=7, T=29, NZ=np.array([0, 15, 14, 11]))
