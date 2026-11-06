"""Static Toffoli-like counts (CCX, CCZ, CSWAP) per named stage of a circuit.

Stages come from the `@enter:` / `@exit:` INFO_TAGs that `builder.mark(...)`
writes. Streams a .kmx text dump, so it works on the 86M-op dqi circuit.
Counts are static (every gate once), not per-shot means.

    python -I tools/profile_marks.py circuit.kmb [DEPTH]
"""
import collections
import os
import sys

import kickmix as km

path = sys.argv[1]
depth = int(sys.argv[2]) if len(sys.argv) > 2 else 3
kmx = path if path.endswith('.kmx') else path + '.tmp.kmx'
if kmx != path:
    km.Circuit.from_file(path).to_file(kmx, format='kmx')

stack, counts = [], collections.Counter()
try:
    with open(kmx) as f:
        for line in f:
            parts = line.split(None, 1)
            if not parts:
                continue
            op = parts[0]
            if op == 'INFO_TAG' and len(parts) > 1:
                tag = parts[1].strip().strip('"')
                if tag.startswith('@enter:'):
                    stack.append(tag[7:])
                elif tag.startswith('@exit:'):
                    stack.pop()
            elif op in ('CCX', 'CCZ', 'CSWAP'):
                n = sum(t.startswith('q') for t in line.split()[1:]) // 3
                counts['/'.join(stack[:depth]) or '(top)'] += max(n, 1)
finally:
    if kmx != path:
        os.remove(kmx)

total = sum(counts.values())
for k, v in sorted(counts.items(), key=lambda kv: -kv[1]):
    print(f'{v:>10,} {100 * v / total:5.1f}%  {k}')
print(f'{total:>10,} total')
