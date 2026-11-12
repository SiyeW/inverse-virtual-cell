"""Render every record in the exploratory one-gene Mixscale projections."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, norm

from render_share_ready_figures import configure

TARGETS = {'TYMSOS': ('ENSG00000176912', 'LINC02315', [18, 223, 108, 275, 57], 9),
           'UBA1': ('ENSG00000130985', 'ZFP42', [16, 86, 99, 123, 214], 31)}
CONTROLS = [('non_targeting_control', 'NTC', 19342),
            ('no_guide_detected', 'No sgRNA detected', 9094),
            ('intergenic_control', 'Intergenic', 195)]
NTC = '#AEDBB5'
TARGET = '#7A6352'
INDEL = '#D73936'


def load_case(root, labels, target):
    gene, marker, sizes, support_count = TARGETS[target]
    path = root / f'{target}_minimum1_scores.csv.gz'
    frame = pd.read_csv(path, keep_default_na=False)
    assert frame.cell_barcode.is_unique and np.isfinite(frame.score).all()
    assert frame.minimum_DE_genes.eq(1).all()
    variants = labels['company']['variants']
    supported = set()
    for event in variants:
        if event['target_gene_id'] == gene and event['kind'] in ('INS', 'DEL') and event['guide_alt'] > 0:
            supported.update((event['guide'], cell['cell_barcode']) for cell in event['cell_support'] if cell['alt'] > 0)
    assert len(supported) == support_count
    frame['indel_support'] = [(guide, cell) in supported for guide, cell in zip(frame.guide_id, frame.cell_barcode)]
    assert int(frame.indel_support.sum()) == support_count
    rows = []
    for role, label, size in CONTROLS:
        part = frame.loc[frame.role.eq(role)].copy()
        assert len(part) == size
        rows.append((label, part))
    ntc = rows[0][1].score.to_numpy()
    np.testing.assert_allclose([ntc.mean(), ntc.std(ddof=1)], [0, 1], atol=1e-8)
    for index, size in enumerate(sizes, 1):
        guide = f'{target}-{index}'
        part = frame.loc[frame.guide_id.eq(guide)].copy()
        assert len(part) == size and part.analysis_group.eq(gene).all()
        rows.append((guide, part))
    assert sum(len(part) for _, part in rows) == len(frame)
    return marker, rows


def summarize(target, rows):
    reference = float(np.quantile(rows[0][1].score, .95))
    summary = []
    for name, frame in rows:
        above = frame.score.gt(reference)
        red = frame.indel_support
        summary.append({'target': target, 'group': name, 'records': len(frame),
            'median_score': float(frame.score.median()), 'mean_score': float(frame.score.mean()),
            'NTC_95th_percentile': reference, 'above_NTC_95th': int(above.sum()),
            'percent_above_NTC_95th': float(above.mean() * 100),
            'RNA_indel_support': int(red.sum()), 'indel_above_NTC_95th': int((red & above).sum())})
    return summary


def draw(target, marker, rows, summary, xlim, tail_limit, path):
    fig, (ax, tail) = plt.subplots(1, 2, figsize=(10.1, 5.8), sharey=True,
        gridspec_kw={'width_ratios': [2.7, 1]})
    fig.subplots_adjust(left=.265, right=.965, top=.795, bottom=.145, wspace=.24)
    ax.set_yticks(range(8), [f'{name} (n = {len(frame):,})' for name, frame in rows])
    ax.set_ylim(7.55, -.6)
    ax.set_xlim(xlim)
    tail.set_xlim(0, tail_limit)
    ax.axvline(0, color='#DDDDDD', lw=.6, zorder=0)
    reference = summary[0]['NTC_95th_percentile']
    ax.axvline(reference, color='#888888', lw=.8, ls='--', zorder=0)
    grid = np.linspace(*xlim, 900)
    total, red_total = 0, 0
    for index, ((name, frame), record) in enumerate(zip(rows, summary)):
        values = frame.score.to_numpy(float)
        deviation = values.std(ddof=1)
        density = gaussian_kde(values, bw_method=.25 / deviation)(grid) if deviation > 0 else norm.pdf(grid, values[0], .25)
        width = .30 * density / density.max()
        support = np.flatnonzero(density > .002 * density.max())
        extent = slice(support[0], support[-1] + 1)
        color = NTC if index == 0 else ('#B5BBC0' if index < 3 else TARGET)
        fill = '#E4F0E6' if index == 0 else ('#F0F1F2' if index < 3 else '#EEE8E3')
        ax.fill_between(grid[extent], index - width[extent], index + width[extent], color=fill, linewidth=0)
        ax.plot(grid[extent], index - width[extent], color=color, lw=.6)
        ax.plot(grid[extent], index + width[extent], color=color, lw=.6)
        jitter = np.random.default_rng(1130 + index).uniform(-.085, .085, len(frame))
        if index < 3:
            points = ax.scatter(values, index + jitter, s=.7, color=color, edgecolors='none', zorder=2)
        else:
            points = ax.scatter(values, index + jitter, s=7, color=color, marker='x', linewidths=.5, zorder=3)
        assert len(points.get_offsets()) == len(frame)
        red = frame.indel_support.to_numpy(bool)
        highlighted = ax.scatter(values[red], index + jitter[red], color=INDEL, s=18, edgecolors='none', zorder=5)
        assert len(highlighted.get_offsets()) == int(red.sum())
        total += len(frame)
        red_total += int(red.sum())
        percent = record['percent_above_NTC_95th']
        tail.barh(index, percent, height=.43, color=color)
        text = f"{record['above_NTC_95th']}/{len(frame):,}\n({percent:.1f}%)"
        if percent > .8 * tail_limit:
            tail.text(percent - .025 * tail_limit, index, text, ha='right', va='center', fontsize=7.7, color='white')
        else:
            tail.text(percent + .035 * tail_limit, index, text, ha='left', va='center', fontsize=7.7)
    ax.axhline(2.5, color='#D7D7D7', lw=.6)
    tail.axhline(2.5, color='#D7D7D7', lw=.6)
    ax.set_xlabel('Single-gene score (NTC mean = 0, SD = 1)')
    tail.set_xlabel('Above NTC 95th percentile\n(% of records)')
    tail.xaxis.set_major_locator(MaxNLocator(4))
    ax.xaxis.set_major_locator(MaxNLocator(6))
    tail.tick_params(axis='y', left=False, labelleft=False)
    for item in [ax, tail]:
        item.spines[['top', 'right', 'left']].set_visible(False)
        item.tick_params(axis='y', length=0)
    fig.suptitle(f'{target} · All chips · minimum DE genes = 1', fontsize=12, y=.975)
    fig.text(.265, .919, f'Exploratory projection using one scoring gene: {marker}', fontsize=9.5)
    handles = [Line2D([], [], color=TARGET, marker='x', ls='none', markersize=4, label='Targeting records'),
               Line2D([], [], color=INDEL, marker='o', ls='none', markersize=4, label='RNA indel support'),
               Line2D([], [], color='#888888', ls='--', lw=.8, label='NTC 95th percentile')]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.258, .895), ncol=3,
               frameon=False, fontsize=8.2, columnspacing=1.3, handletextpad=.5)
    fig.text(.265, .045, 'The NTC percentile is a display reference; no Strong or KO calls are assigned.', fontsize=8.5, color='#555555')
    fig.savefig(path, dpi=300)
    axes = [list(item.get_position().bounds) for item in [ax, tail]]
    plt.close(fig)
    return {'target': target, 'marker': marker, 'all_score_records_shown': total,
            'RNA_indel_support_shown': red_total, 'axes': axes,
            'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--labels', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    labels = json.loads(args.labels.read_bytes())
    cases = {target: load_case(args.input, labels, target) for target in TARGETS}
    summaries = {target: summarize(target, rows) for target, (_, rows) in cases.items()}
    values = np.concatenate([frame.score.to_numpy() for _, rows in cases.values() for _, frame in rows])
    margin = .03 * np.ptp(values)
    xlim = (float(values.min() - margin), float(values.max() + margin))
    maximum = max(row['percent_above_NTC_95th'] for summary in summaries.values() for row in summary)
    tail_limit = min(100, max(20, np.ceil(maximum * 1.45 / 5) * 5))
    configure(Path('C:/Windows/Fonts/arial.ttf'))
    args.output.mkdir(parents=True, exist_ok=False)
    receipts = [draw(target, marker, rows, summaries[target], xlim, tail_limit,
                     args.output / f'{target}_minimum1_scores.png') for target, (marker, rows) in cases.items()]
    assert receipts[0]['axes'] == receipts[1]['axes']
    pd.DataFrame([row for summary in summaries.values() for row in summary]).to_csv(args.output / 'group_score_summary.csv', index=False)
    receipt = {'cases': receipts, 'score_limits': xlim, 'percent_limits': [0, tail_limit],
        'source_sha256': {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [args.input / 'TYMSOS_minimum1_scores.csv.gz', args.input / 'UBA1_minimum1_scores.csv.gz', args.labels]}}
    (args.output / 'render_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps({'cases': receipts, 'summary': summaries}, indent=2))


if __name__ == '__main__':
    main()
