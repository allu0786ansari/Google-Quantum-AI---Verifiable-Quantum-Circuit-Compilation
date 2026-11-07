"""Build the 893-qubit dqi decoder from ahoirg's authenticated source bundle, optionally
with a builder proxy that (a) tags frees and/or (b) rewrites gates on the fly.

    python build_tagged.py BUNDLE WORK OUT.kmb [--tag-frees]

Mirrors build_dqi.py (same imports, same transform.install, same register metadata).
"""
import argparse, hashlib, importlib.util, json, sys, time
from pathlib import Path, PurePosixPath

ap = argparse.ArgumentParser()
ap.add_argument('bundle', type=Path); ap.add_argument('work', type=Path); ap.add_argument('out')
ap.add_argument('--tag-frees', action='store_true')
args = ap.parse_args()

import kickmix as km
t0 = time.time()
bundle = json.loads(args.bundle.read_text())
work = args.work.resolve()
if not work.exists():
    for name, pin in bundle['files'].items():
        rel = PurePosixPath(name)
        assert not rel.is_absolute() and '..' not in rel.parts
        data = pin['text'].encode()
        assert hashlib.sha256(data).hexdigest() == pin['sha256']
        (work / name).parent.mkdir(parents=True, exist_ok=True)
        (work / name).write_bytes(data)
source = work / 'blocks/dqi/pro-q18-900/source'
helpers = work / 'portable_builder'
sys.path.insert(0, str(helpers))
sys.modules['circ'] = km
spec = importlib.util.spec_from_file_location('carry_decoder', source / 'decode_best.py')
module = importlib.util.module_from_spec(spec); sys.modules[spec.name] = module; spec.loader.exec_module(module)
sys.path.insert(0, str(helpers))
from transform import install
manager, transformation = install(module.decode, source)
print('imported', round(time.time() - t0), 's', flush=True)


class Proxy:
    """Delegates to a real CircuitBuilder; tags frees of qubits."""
    def __init__(self, b):
        object.__setattr__(self, '_b', b)

    def __getattr__(self, name):
        return getattr(self._b, name)

    def __setattr__(self, name, value):
        setattr(self._b, name, value)

    def free(self, arg):
        b = self._b
        items = arg if isinstance(arg, (list, tuple)) else ([arg] if not hasattr(arg, '__len__') else list(arg))
        qs = [x for x in items if isinstance(x, km.q)]
        if qs:
            b.info_tag('@free ' + ' '.join(str(q) for q in qs))
        return b.free(arg)


real = km.CircuitBuilder()
builder = Proxy(real) if args.tag_frees else real
syndromes = [real.alloc_qubits(12) for _ in range(70)]
masks = [real.alloc_bits(12) for _ in range(4095)]
module.decode(builder, syndromes, masks)
meta = []
for i in range(70):
    meta.append(f'REGISTER r{i} "S{i}"'); meta.extend(f'APPEND_TO_REGISTER q{12 * i + k} r{i}' for k in range(12))
for j in range(4095):
    meta.append(f'REGISTER r{70 + j} "E{j}"'); meta.extend(f'APPEND_TO_REGISTER b{12 * j + k} r{70 + j}' for k in range(12))
circuit = km.Circuit('\n'.join(meta)) + real.finish_circuit()
circuit.to_file(args.out, format='kmb' if args.out.endswith('.kmb') else 'kmx')
print(json.dumps(dict(out=args.out, Q=circuit.num_qubits, ops=len(circuit),
                      static_toffoli=circuit.operation_statistics().num_toffoli_like,
                      seconds=round(time.time() - t0))), flush=True)
