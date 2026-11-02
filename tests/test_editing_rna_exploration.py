import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from explore_editing_rna_cells import cut_position, family_state, covariance_pca
from analyze_guided_rna_alleles import Read, raw_events, covers_reference


def read(seq, cigar, start=100):
    sam = '\t'.join(['r', '0', 'chr1', str(start + 1), '60', cigar, '*', '0', '0',
                     seq, 'I' * len(seq), 'CB:Z:c1', 'UB:Z:u1', 'GX:Z:g1'])
    return Read(sam, 'G1')


def test_checked_plus_and_minus_canonical_cut_sites():
    spacer = 'A' * 20
    assert cut_position(spacer, spacer, 100, lambda a, b: 'AGG')['cut0'] == 117
    result = cut_position(spacer, 'T' * 20, 100, lambda a, b: 'CCA')
    assert result['strand'] == '-' and result['cut0'] == 103
    assert cut_position(spacer, spacer, 100, lambda a, b: 'AAA')['cut0'] is None


def test_splice_skip_is_not_a_deletion_or_reference_coverage():
    r = read('A' * 10, '5M100N5M')
    region = {'guides': [{'guide': 'g', 'start0': 100, 'end0': 215}]}
    assert raw_events(r, region, lambda c, a, b: 'A' * (b - a)) == []
    assert covers_reference(r, 102, 210) is None


def test_pcr_disagreement_is_retained_as_ambiguous():
    r = read('A' * 20, '20M')
    variant = ('chr1', 105, 'A', 'T')
    assert family_state([(r, {variant}), (r, set())], variant) == 'ambiguous'
    assert family_state([(r, {variant}), (r, {variant})], variant) == 'alt'
    assert family_state([(r, set())], variant) == 'ref'


def test_uncovered_variant_is_unknown_not_reference():
    r = read('A' * 5, '5M')
    assert family_state([(r, set())], ('chr1', 120, 'A', 'T')) is None


def test_all_record_covariance_pca_matches_direct_svd_geometry():
    rng = np.random.default_rng(71)
    ntc = rng.normal(size=(13, 5)); query = rng.normal(size=(8, 5))
    bp, qp, fraction = covariance_pca(ntc, query, components=5)
    x = np.concatenate([ntc, query]); x -= x.mean(axis=0)
    scores = np.concatenate([bp, qp])
    assert np.allclose(scores @ scores.T, x @ x.T)
    assert np.allclose(fraction.sum(), 1)
