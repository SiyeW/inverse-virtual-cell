"""Create a small paired FASTQ prefix for engineering checks, never inference."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import subprocess

from manage_core_study import input_inventory


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def read_record(stream):
    record = [stream.readline() for _ in range(4)]
    if not record[0]:
        raise ValueError('Source ended before requested prefix')
    if any(not x.endswith(b'\n') for x in record) or not record[0].startswith(b'@') or not record[2].startswith(b'+'):
        raise ValueError('Malformed FASTQ record')
    if len(record[1].rstrip(b'\r\n')) != len(record[3].rstrip(b'\r\n')):
        raise ValueError('FASTQ sequence/quality lengths differ')
    return record


def read_id(record):
    value = record[0].split()[0]
    return value[:-2] if value.endswith((b'/1', b'/2')) else value


def copy_pair(left, right, out_left, out_right, pairs):
    if type(pairs) is not int or not 1 <= pairs <= 200000:
        raise ValueError('Engineering prefix requires 1..200000 read pairs')
    sources, outputs = [Path(left), Path(right)], [Path(out_left), Path(out_right)]
    if len({p.resolve() for p in sources + outputs}) != 4 or any(p.exists() for p in outputs):
        raise ValueError('Source/output overlap or prior output')
    before = [p.stat() for p in sources]
    h = [hashlib.sha256(), hashlib.sha256()]
    unpacked = [0, 0]
    with gzip.open(sources[0], 'rb') as a, gzip.open(sources[1], 'rb') as b, outputs[0].open('xb') as oa, outputs[1].open('xb') as ob:
        with gzip.GzipFile(filename='', mode='wb', fileobj=oa, mtime=0) as za, gzip.GzipFile(filename='', mode='wb', fileobj=ob, mtime=0) as zb:
            for _ in range(pairs):
                ra, rb = read_record(a), read_record(b)
                if read_id(ra) != read_id(rb):
                    raise ValueError('R1/R2 read IDs differ')
                for i, (record, stream) in enumerate([(ra, za), (rb, zb)]):
                    raw = b''.join(record);h[i].update(raw);unpacked[i] += len(raw)
                    if sum(unpacked) > 256 * 1024**2:
                        raise ValueError('Engineering prefix byte budget exceeded')
                    stream.write(raw)
    results = []
    for i, p in enumerate(sources):
        after = p.stat()
        if (before[i].st_size, before[i].st_mtime_ns, before[i].st_dev) != (after.st_size, after.st_mtime_ns, after.st_dev):
            raise ValueError('Source stable identity changed')
        actual = hashlib.sha256()
        with gzip.open(outputs[i], 'rb') as stream:
            for block in iter(lambda: stream.read(4 * 1024**2), b''):
                actual.update(block)
        if actual.hexdigest() != h[i].hexdigest():
            raise ValueError('Written prefix differs from source records')
        results.append(dict(source=str(p), output=str(outputs[i]), pairs=pairs,
            source_prefix_uncompressed_sha256=h[i].hexdigest(), output_sha256=sha256(outputs[i]),
            output_bytes=outputs[i].stat().st_size, uncompressed_bytes=unpacked[i]))
    return results


def prepare(config, chip_id, pairs, output, allowed_parent):
    output, allowed_parent = Path(output).resolve(), Path(allowed_parent).resolve()
    if output == allowed_parent or not output.is_relative_to(allowed_parent):
        raise ValueError('Output must be a new child of the allowed parent')
    if any(output == Path(p).resolve() or output.is_relative_to(Path(p).resolve()) for p in config['protected_roots']):
        raise ValueError('Output overlaps protected input')
    selected = [c for c in config['chips'] if c['chip_id'] == chip_id]
    if len(selected) != 1:
        raise ValueError('Exactly one known chip is required')
    if len(selected[0]['libraries']) != 2:
        raise ValueError('This bounded prefix entry requires one GEX and one CRISPR library')
    if {lib['feature_type'] for lib in selected[0]['libraries']} != {'Gene Expression', 'CRISPR Guide Capture'}:
        raise ValueError('One GEX and one CRISPR library required')
    subset = dict(config, chips=selected)
    inventory = input_inventory(subset)  # old full receipts + current metadata, no bulk rehash
    by_path = {r['path']: r for r in inventory}
    output.mkdir(exist_ok=False)
    records = []
    for i, library in enumerate(selected[0]['libraries']):
        if library['feature_type'] not in ['Gene Expression', 'CRISPR Guide Capture']:
            raise ValueError('Unsupported engineering library')
        keys = library['receipt_keys']
        r1, r2 = [[Path(library['fastq_dir']) / Path(k).name for k in keys if f'_R{read}_' in k] for read in (1, 2)]
        if len(r1) != 1 or len(r2) != 1:
            raise ValueError('One unambiguous R1/R2 pair required per library')
        dest = output / f'library_{i}';dest.mkdir()
        copied = copy_pair(r1[0], r2[0], dest/r1[0].name, dest/r2[0].name, pairs)
        for row in copied:
            row['prior_full_source_receipt'] = by_path[row['source']]
        records.append(dict(fastq_id=library['fastq_id'], feature_type=library['feature_type'], fastq_dir=str(dest), files=copied))
    receipt = dict(scientific_result=False, sampling='first_read_pairs_not_representative',
                   chip_id=chip_id, pairs_per_library=pairs, libraries=records)
    (output/'prefix_receipt.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf8')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['config','chip-id','output','allowed-parent','expected-commit']:
        parser.add_argument('--'+name, required=True)
    parser.add_argument('--pairs', type=int, default=100000)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    commit = subprocess.check_output(['git','rev-parse','HEAD'], cwd=repo, text=True).strip()
    if commit != args.expected_commit or subprocess.check_output(['git','status','--porcelain'], cwd=repo, text=True).strip():
        raise ValueError('Use a fixed clean pushed execution revision')
    result = prepare(json.loads(Path(args.config).read_text()), args.chip_id, args.pairs, args.output, args.allowed_parent)
    print(json.dumps(dict(chip_id=result['chip_id'], pairs_per_library=result['pairs_per_library'], scientific_result=False)))


if __name__ == '__main__':main()
