"""Scoring harness: exact scores, ratios and cost breakdowns for all three tasks.

    export VERIFIERS=/path/to/verifiers
    python -I tools/harness.py invsqrt circuit.kmb [--max-x 65536] [--full]
    python -I tools/harness.py dqi circuit.kmb [--shots 1024]
    python -I tools/harness.py chem submission_dir/
    python -I tools/harness.py chem-model R B C BROT --angles 8932 --lam 47.03 --J 5
    python -I tools/harness.py profile circuit.kmb [--depth 2]
    python -I tools/harness.py overall 21.098 65.106 10.698

Scores use the official verifier code. Quick modes check fewer inputs/shots; the
score they print is an estimate (invsqrt: rotation precision depends on the worst
margin over all inputs; dqi: assumes 0 failures over the 102,400 leaderboard shots).
"""
import argparse
import collections
import json
import math
import os
import sys
import tempfile
import time
import types

VERIFIERS = os.environ.get('VERIFIERS', 'verifiers')
sys.path.insert(0, VERIFIERS)

BASELINE = {'invsqrt': 2076774.23, 'dqi': 24439397322564.6, 'dfthcsos': 3355599588.3155375}
FRONTIER = {'invsqrt': 21.098, 'dqi': 65.106, 'dfthcsos': 10.698}
DQI_R_ZERO_FAILURES = 1.19126   # Clopper-Pearson factor at 0 failures / 102,400 shots, delta 2^-128


def overall(ratios):
    return math.prod(ratios) ** (1 / len(ratios))


def what_if(task, ratio):
    r = dict(FRONTIER, **{task: ratio})
    return overall(r.values())


def score_invsqrt(path, max_x=65536, full=False):
    import kickmix as km
    from invsqrt import verify_submission_invsqrt as V
    c = km.Circuit.from_file(path)
    spec = (V.ProblemSpecInvSqrt.for_n_p(10) if full else
            V.ProblemSpecInvSqrt(width=22, max_x=max_x, beta=V.Fraction(1), precision_bits=24,
                                 shots_per_x=4, system_qubits=0))
    r = V.verify(c, spec, batch_size=256, max_workers=os.cpu_count())
    keys = ('status', 'score', 'num_qubits', 'toffoli_count', 't_count', 'rotation_count',
            'rotation_eps', 'toffoli_equivalent', 'worst_error_over_tolerance')
    out = {k: r[k] for k in keys if k in r}
    out['ratio'] = BASELINE['invsqrt'] / r['score']
    return out


def _load_dqi_verifier():
    src_path = os.path.join(VERIFIERS, 'dqi', 'verify_submission_dqi.py')
    src = open(src_path).read().replace('circuit.reaction_depth()', '-1')   # report only; OOMs
    V = types.ModuleType('vdqi')
    V.__file__ = src_path
    sys.modules['vdqi'] = V
    exec(compile(src, src_path, 'exec'), V.__dict__)
    return V


def score_dqi(path, shots=1024):
    import kickmix as km
    V = _load_dqi_verifier()
    c = km.Circuit.from_file(path)
    r = V.verify(c, V.ProblemSpecDQI.benchmark_instance(num_shots=shots), max_workers=1, log2_delta=-1.0)
    ops = r.pop('shot_operation_counts', {})
    q = r['num_qubits']
    t = r['toffoli_count']
    est = DQI_R_ZERO_FAILURES * q * q * t
    return dict(status=r.get('status'), failures=r.get('num_failures', r.get('failures')), num_qubits=q,
                mean_toffolis=t, score_if_0_failures=est, ratio=BASELINE['dqi'] / est,
                ops_per_shot={k: v / shots for k, v in ops.items()})


def score_chem(folder):
    from dfthcsos.verify import verify
    r = verify(str(folder))
    res, m = r['resources'], r['metrics']
    return dict(status=r['status'], score=r['score'], ratio=BASELINE['dfthcsos'] / r['score'],
                lambda_eff=r['lambda_eff'], toffolis=r['toffoli_count'], qubits=r['num_qubits'],
                brot=r['brot'], J=r['batch_size'], K=r['num_circuits'],
                angles_toffoli=res['angles_toffoli'], fixed_toffoli=res['heuristic_toffoli'],
                persistent=res['persistent_qubits'], temporary=res['temporary_qubits'],
                error_1_mHa=1e3 * m['error_1'], error_2_mHa=1e3 * m['error_2'],
                Lambda=m['Lambda'], E_gap=m['E_gap'])


def chem_model(R, B, C, brot, angles_toffoli, lam, J, N=76):
    """lambda * T * Q from the verifier's own heuristic cost model plus the angle loaders."""
    from dfthcsos.verify import heuristic_costs
    h = heuristic_costs(N, R, B, C, brot)
    T = angles_toffoli + h['toffoli']
    Q = h['persistent_qubits'] + J * brot + h['temporary_qubits']
    score = lam * T * Q
    return dict(R=R, B=B, C=C, brot=brot, J=J, vectors=N + R * B, fixed_toffoli=h['toffoli'],
                T=T, Q=Q, inner_temp=h['inner_state_prep_temporary_qubits'],
                outer_temp=h['outer_state_prep_temporary_qubits'], score=score,
                ratio=BASELINE['dfthcsos'] / score)


def profile(path, depth=2):
    """Toffoli-like gates per stage tag; conditioned ones (fire on ~half the shots) apart."""
    import kickmix as km
    with tempfile.TemporaryDirectory() as tmp:
        kmx = os.path.join(tmp, 'c.kmx')
        km.Circuit.from_file(path).to_file(kmx, format='kmx')
        stack, counts = [], collections.defaultdict(collections.Counter)
        with open(kmx) as f:
            for line in f:
                parts = line.split(None, 1)
                if not parts:
                    continue
                op = parts[0]
                if op == 'INFO_TAG':
                    tag = parts[1].strip().strip('"')
                    if tag.startswith('@enter:'):
                        stack.append(tag[7:])
                    elif tag.startswith('@exit:'):
                        stack.pop()
                    continue
                if op in ('CCX', 'CCZ', 'CSWAP', 'HMR'):
                    n = sum(t.startswith('q') for t in line.split()[1:])
                    n = max(1, n // 3) if op != 'HMR' else max(1, n)
                    key = op + ('?' if ' if ' in line else '')
                    counts['/'.join(stack[:depth]) or '(top)'][key] += n
    rows = []
    for stage, c in counts.items():
        tof = sum(v for k, v in c.items() if k[:3] in ('CCX', 'CCZ', 'CSW'))
        rows.append((tof, stage, dict(c)))
    return sorted(rows, reverse=True)


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest='cmd', required=True)
    a = sub.add_parser('invsqrt'); a.add_argument('path'); a.add_argument('--max-x', type=int, default=65536)
    a.add_argument('--full', action='store_true')
    a = sub.add_parser('dqi'); a.add_argument('path'); a.add_argument('--shots', type=int, default=1024)
    a = sub.add_parser('chem'); a.add_argument('path')
    a = sub.add_parser('chem-model')
    for k in ('R', 'B', 'C', 'brot'):
        a.add_argument(k, type=int)
    a.add_argument('--angles', type=float, default=8932); a.add_argument('--lam', type=float, default=47.0256)
    a.add_argument('--J', type=int, default=5)
    a = sub.add_parser('profile'); a.add_argument('path'); a.add_argument('--depth', type=int, default=2)
    a = sub.add_parser('overall'); a.add_argument('ratios', type=float, nargs=3)
    args = p.parse_args()
    t0 = time.time()
    if args.cmd == 'invsqrt':
        out = score_invsqrt(args.path, args.max_x, args.full)
        out['overall_with_frontier_rest'] = what_if('invsqrt', out['ratio'])
    elif args.cmd == 'dqi':
        out = score_dqi(args.path, args.shots)
        out['overall_with_frontier_rest'] = what_if('dqi', out['ratio'])
    elif args.cmd == 'chem':
        out = score_chem(args.path)
        out['overall_with_frontier_rest'] = what_if('dfthcsos', out['ratio'])
    elif args.cmd == 'chem-model':
        out = chem_model(args.R, args.B, args.C, args.brot, args.angles, args.lam, args.J)
        out['overall_with_frontier_rest'] = what_if('dfthcsos', out['ratio'])
    elif args.cmd == 'profile':
        for tof, stage, c in profile(args.path, args.depth)[:40]:
            print(f'{tof:>9,}  {stage[:60]:60s} {c}')
        return
    else:
        out = {'overall': overall(args.ratios)}
    print(json.dumps(out, indent=1, default=float))
    print(f'({time.time() - t0:.1f}s)', file=sys.stderr)


if __name__ == '__main__':
    main()
