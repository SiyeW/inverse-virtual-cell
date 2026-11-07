"""Highlight the three observed indels of one guide on accepted RNA coordinates."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from render_share_ready_figures import configure, save, assert_drawn_points


def read_json(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('coordinates', 'coordinate-receipt', 'sequences', 'coding', 'output'):
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--guide', required=True)
    parser.add_argument('--commit', required=True)
    args = parser.parse_args()
    receipt = read_json(args.coordinate_receipt)
    assert sha(args.coordinates) == receipt['files'][args.coordinates.name]['sha256']
    source = read_json(args.sequences)['company']
    coding = read_json(args.coding)
    changes = [v for v in source['variants'] if v['guide'] == args.guide and v['kind'] in ('DEL', 'INS')]
    changes.sort(key=lambda v: (v['kind'] != 'DEL', v['length_change']))
    assert [(v['kind'], v['length_change']) for v in changes] == [('DEL', -9), ('INS', 1), ('INS', 8)]
    target = args.guide.rsplit('-', 1)[0]
    with np.load(args.coordinates, allow_pickle=False) as arrays:
        xy = arrays['UMAP'].copy()
        guides = arrays['guide'].astype(str)
        barcodes = arrays['cell_barcode'].astype(str)
        roles = arrays['role'].astype(str)
    assert xy.shape == (len(guides), 2) and np.isfinite(xy).all()
    keys = list(zip(guides, barcodes))
    assert len(set(keys)) == len(keys)
    lookup = {key: i for i, key in enumerate(keys)}
    target_mask = np.array([g.rsplit('-', 1)[0] == target for g in guides])
    ntc_mask = roles == 'non_targeting_control'
    assert target_mask.sum() == 537 and ntc_mask.sum() == 19342
    assert (guides == args.guide).sum() == 184
    background = xy[target_mask | ntc_mask]
    center = (background.min(axis=0) + background.max(axis=0)) / 2
    half_span = np.ptp(background, axis=0).max() * .53
    configure(Path('C:/Windows/Fonts/arial.ttf'))
    fig, axes = plt.subplots(1, 3, figsize=(10.1, 3.9), sharex=True, sharey=True)
    fig.subplots_adjust(left=.065, right=.99, bottom=.19, top=.88, wspace=.12)
    mapping = []
    for ax, variant, n in zip(axes, changes, (2, 16, 1)):
        selected_barcodes = sorted({r['cell_barcode'] for r in variant['cell_support'] if r['alt'] > 0})
        assert len(selected_barcodes) == n
        indices = [lookup[(args.guide, barcode)] for barcode in selected_barcodes]
        highlights = xy[indices]
        event = next(e for e in coding['events'] if all(e[k] == variant[k]
                     for k in ('guide', 'chrom', 'pos0', 'ref', 'alt')))
        assert event['category'] == ('In-frame' if variant['kind'] == 'DEL' else 'Frameshift')
        ax.scatter(xy[ntc_mask, 0], xy[ntc_mask, 1], s=1.2, color='#496E8B', alpha=.10,
                   edgecolors='none', rasterized=True, zorder=1)
        ax.scatter(xy[target_mask, 0], xy[target_mask, 1], s=6, color='#B5BDC3', alpha=.42,
                   edgecolors='none', rasterized=True, zorder=2)
        ax.scatter(highlights[:, 0], highlights[:, 1], s=29, color='#B34B46', alpha=1,
                   edgecolors='white', linewidths=.35, zorder=3)
        assert_drawn_points(ax, np.vstack((xy[ntc_mask], xy[target_mask], highlights)))
        word = 'deletion' if variant['kind'] == 'DEL' else 'insertion'
        ax.set_title(f"{abs(variant['length_change'])} bp {word} (n = {n})")
        ax.set(xlabel='UMAP 1', xlim=(center[0]-half_span, center[0]+half_span),
               ylim=(center[1]-half_span, center[1]+half_span))
        ax.set_aspect('equal', adjustable='box')
        mapping.append(dict(guide=args.guide, chrom=variant['chrom'], pos0=variant['pos0'],
                            ref=variant['ref'], alt=variant['alt'], change_bp=variant['length_change'],
                            category=event['category'], barcode_records=n,
                            RNA_UMI_families=variant['guide_alt'],
                            barcodes=selected_barcodes, coordinates=highlights.tolist()))
    axes[0].set_ylabel('UMAP 2')
    handles = [Line2D([], [], ls='', marker='o', color=color, markersize=size, label=label)
               for color, size, label in [('#496E8B', 3, 'NTC'), ('#B5BDC3', 3, 'All target records'),
                                           ('#B34B46', 4, 'Supports indicated change')]]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.51, .005),
               ncol=3, columnspacing=1.6, handletextpad=.45)
    args.output.mkdir(parents=True, exist_ok=False)
    save(fig, args.output, 'three_indel_RNA_states')
    result = dict(code_commit=args.commit, target=target, guide=args.guide, temporary_display=True,
                  representation='Existing whole-cohort RNA UMAP; not the target-specific PCA in the old HTML',
                  RNA_coordinates_refitted=False, all_NTC_records=19342, all_target_records=537,
                  same_background_and_axes=True, all_points_verified_against_source=True,
                  raw_BAM_or_count_matrix_read=False, changes=mapping,
                  inputs_sha256={str(p): sha(p) for p in [args.coordinates, args.coordinate_receipt,
                                  args.sequences, args.coding]},
                  figures_sha256={p.name: sha(p) for p in args.output.iterdir()})
    (args.output / 'render_receipt.json').write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'highlighted_records': [r['barcode_records'] for r in mapping],
                      'all_target_records': 537, 'all_NTC_records': 19342, 'coordinates_unchanged': True}))


if __name__ == '__main__':
    main()
