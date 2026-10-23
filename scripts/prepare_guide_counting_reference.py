"""Separate guide capture definitions from optional target-expression annotations.

The immutable source remains authoritative for downstream target/control labels.
This tool neither relabels controls nor changes guide detection fields.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess

CAPTURE = ['id', 'name', 'read', 'pattern', 'sequence', 'feature_type']
TARGET = ['target_gene_id', 'target_gene_name']


def digest(data):
    return hashlib.sha256(data).hexdigest()


def encode(rows, fields):
    text = io.StringIO(newline='')
    writer = csv.DictWriter(text, fieldnames=fields, lineterminator='\n')
    writer.writeheader()
    writer.writerows([{k: row[k] for k in fields} for row in rows])
    return text.getvalue().encode('utf8')


def prepare(source, expected_sha256, output, allowed_parent, code=None):
    if any(not Path(p).is_absolute() for p in (source, output, allowed_parent)):
        raise ValueError('Absolute source and output paths required')
    source, output, allowed_parent = (Path(p).resolve() for p in (source, output, allowed_parent))
    if not output.is_relative_to(allowed_parent) or output == allowed_parent:
        raise ValueError('Output must be a new child of the allowed parent')
    if source.is_relative_to(output) or output.is_relative_to(source.parent):
        raise ValueError('Output overlaps source delivery directory')
    if not re.fullmatch('[0-9a-f]{64}', expected_sha256) or source.stat().st_size > 8 * 1024**2:
        raise ValueError('Expected SHA256 and bounded source required')
    before = source.stat()
    raw = source.read_bytes()
    if digest(raw) != expected_sha256:
        raise ValueError('Source SHA256 differs')
    reader = csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))
    if reader.fieldnames is None or len(reader.fieldnames) != len(set(reader.fieldnames)) or set(reader.fieldnames) != set(CAPTURE + TARGET):
        raise ValueError('Expected exactly six capture and two target columns')
    rows = list(reader)
    if not rows or any(set(row) != set(reader.fieldnames) or any(row[k] is None for k in reader.fieldnames) for row in rows):
        raise ValueError('Malformed or empty feature reference')
    if len({r['id'] for r in rows}) != len(rows):
        raise ValueError('Duplicate guide IDs')
    for row in rows:
        if any(not row[k] for k in CAPTURE) or row['feature_type'] != 'CRISPR Guide Capture':
            raise ValueError('Complete CRISPR capture fields required')
        if row['read'] not in {'R1', 'R2'} or '(BC)' not in row['pattern'] or not re.fullmatch('[ACGTN]+', row['sequence']):
            raise ValueError('Invalid read, capture pattern or sequence')
    payloads = {
        'source_feature_reference.csv': raw,
        'guide_counting_features.csv': encode(rows, CAPTURE),
        'authoritative_target_map.csv': encode(rows, ['id', 'name'] + TARGET),
    }
    # The derived capture table must retain every guide and every capture field,
    # in the same order. Missing target annotations never become Non-Targeting.
    reread = list(csv.DictReader(io.StringIO(payloads['guide_counting_features.csv'].decode())))
    if reread != [{k:r[k] for k in CAPTURE} for r in rows]:
        raise ValueError('Capture fields changed during serialization')
    after = source.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns) or source.read_bytes() != raw:
        raise ValueError('Source changed while preparing')
    output.mkdir(parents=False, exist_ok=False)
    for name, data in payloads.items():
        (output / name).write_bytes(data)
    receipt = dict(schema_version=1, purpose='guide_counting_without_optional_target_expression_annotations',
        source=str(source), source_sha256=expected_sha256, guides=len(rows),
        capture_fields_unchanged=True, target_and_control_relabeling=False,
        downstream_target_authority='source_feature_reference.csv and authoritative_target_map.csv; retain separately verified control catalogue',
        upstream_target_expression_annotation_enabled=False, scientific_result=False,
        code=code, files={k:dict(bytes=len(v), sha256=digest(v)) for k,v in payloads.items()})
    (output / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf8')
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('source', 'source-sha256', 'output', 'allowed-parent', 'expected-commit'):
        parser.add_argument('--' + key, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip()
    dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=repo, text=True).strip()
    if commit != args.expected_commit or dirty:
        raise ValueError('Use the pushed, clean, fixed execution revision')
    print(json.dumps(prepare(args.source, args.source_sha256, args.output, args.allowed_parent, dict(commit=commit))))


if __name__ == '__main__':
    main()
