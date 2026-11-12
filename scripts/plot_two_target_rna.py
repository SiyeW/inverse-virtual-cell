"""Plot complete per-guide Mixscale distributions and SCEPTRE RNA counts."""
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

from plot_indel_pca_by_target import IND_WITH_NTC, SAME
from render_share_ready_figures import configure


EXPECTED = {'TYMSOS': ('ENSG00000176912', 681, 9), 'UBA1': ('ENSG00000130985', 538, 31)}
CONTROLS = [('non_targeting_control', 'NTC', 19342),
            ('no_guide_detected', 'No sgRNA detected', 9094),
            ('intergenic_control', 'Intergenic', 195)]


def load_inputs(root, labels_path, sceptre_stage='sceptre'):
    labels = json.loads(labels_path.read_bytes())['company']['variants']
    counts = pd.read_csv(root / sceptre_stage / 'guide_DE_counts.csv', keep_default_na=False)
    if not counts.guide_id.is_unique or len(counts) != 10:
        raise ValueError('Expected every guide in the two-target pilot')
    cases = {}
    for target, (gene_id, expected_n, expected_support) in EXPECTED.items():
        frame = pd.read_csv(root / f'mixscale/{target}_scores.csv.gz', keep_default_na=False)
        frame['score'] = pd.to_numeric(frame.score, errors='coerce')
        query = frame.analysis_group.eq(gene_id)
        guides = sorted(frame.loc[query, 'guide_id'].unique())
        assert guides == [target + '-' + str(i) for i in range(1, 6)]
        assert int(query.sum()) == expected_n and frame.cell_barcode.is_unique
        supported = set()
        for event in labels:
            if event['target_gene_id'] == gene_id and event['kind'] in ('INS', 'DEL') and event['guide_alt'] > 0:
                supported.update((event['guide'], cell['cell_barcode']) for cell in event['cell_support'] if cell['alt'] > 0)
        assert len(supported) == expected_support
        present = set(zip(frame.loc[query, 'guide_id'], frame.loc[query, 'cell_barcode']))
        assert supported <= present
        frame['indel_support'] = [(guide, cell) in supported for guide, cell in zip(frame.guide_id, frame.cell_barcode)]
        rows = []
        for role, label, number in CONTROLS:
            part = frame.loc[frame.role.eq(role)]
            assert len(part) == number
            rows.append((label, part))
        for guide in guides:
            part = frame.loc[frame.guide_id.eq(guide)]
            reference = counts.loc[counts.guide_id.eq(guide)].iloc[0]
            assert len(part) == int(reference.records)
            rows.append((guide, part))
        finite_ntc = rows[0][1].score.dropna().to_numpy()
        if len(finite_ntc):
            assert len(finite_ntc) == 19342
            np.testing.assert_allclose([finite_ntc.mean(), finite_ntc.std(ddof=1)], [0, 1], atol=1e-6)
        cases[target] = (rows, counts.loc[counts.guide_id.isin(guides)].set_index('guide_id'))
    return cases


def draw_unscored(target, rows, counts, bar_limit, destination, feature_count):
    """Show the observed tests and coverage when the score model is undefined."""
    if feature_count is None or feature_count >= 5:
        raise ValueError('Missing scores require the observed feature-gate result')
    guides = rows[3:]
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    fig.subplots_adjust(left=.53, right=.95, bottom=.17, top=.70)
    y = np.arange(len(guides))
    ax.set_ylim(4.6, -.6)
    ax.set_xlim(0, bar_limit)
    ax.set_yticks([])
    fig.suptitle(f'{target} · All chips ({sum(len(frame) for _, frame in guides):,} records)',
                 fontsize=12, y=.965)
    fig.text(.075, .87, f'Mixscale: {feature_count} eligible RNA; at least 5 required.', fontsize=10)
    fig.text(.075, .815, 'Single-cell response scores could not be estimated.', fontsize=10)
    for x, label in [(.075, 'sgRNA'), (.27, 'Records'), (.385, 'RNA indel support')]:
        fig.text(x, .74, label, fontsize=9, ha='center' if x > .1 else 'left')
    positions = []
    for index, (name, frame) in enumerate(guides):
        point = fig.transFigure.inverted().transform(ax.transData.transform((0, index)))
        positions.append(point[1])
        fig.text(.075, point[1], name, fontsize=10, va='center')
        fig.text(.27, point[1], str(len(frame)), fontsize=10, va='center', ha='center')
        support = int(frame.indel_support.sum())
        fig.text(.385, point[1], str(support), fontsize=10, va='center', ha='center',
                 color=IND_WITH_NTC if support else '#777777')
        item = counts.loc[name]
        if item.assessment != 'evaluated':
            ax.text(.02 * bar_limit, index, 'Not assessed', va='center', fontsize=9)
            continue
        value = int(item.DE_RNAs)
        ax.barh(index, value, height=.48, color=SAME)
        ax.text(value + .025 * bar_limit, index, str(value), va='center', fontsize=10)
    ax.set_xlabel('DE RNAs\n(SCEPTRE; FDR < 0.1)', fontsize=10)
    ax.xaxis.set_major_locator(MaxNLocator(4, integer=True))
    ax.spines[['top', 'right', 'left']].set_visible(False)
    fig.text(.075, .07, f'Reference: NTC ({len(rows[0][1]):,} records). All five sgRNAs retained.',
             fontsize=9, color='#555555')
    fig.savefig(destination, dpi=300)
    bounds = list(ax.get_position().bounds)
    plt.close(fig)
    return {'target': target, 'score_model': 'not_estimable', 'eligible_features': feature_count,
            'required_features': 5, 'finite_score_records': 0, 'plotted_indel_records': 0,
            'indel_support_records_shown': sum(int(frame.indel_support.sum()) for _, frame in guides),
            'rows': [{'label': name, 'records': len(frame), 'scored': 0} for name, frame in rows],
            'axes': [bounds], 'sha256': hashlib.sha256(destination.read_bytes()).hexdigest()}


def draw(target, rows, counts, xlim, bar_limit, destination, feature_count=None):
    if not any(frame.score.notna().any() for _, frame in rows):
        return draw_unscored(target, rows, counts, bar_limit, destination, feature_count)
    fig, (score_ax, count_ax) = plt.subplots(1, 2, figsize=(9.2, 5.0), sharey=True,
        gridspec_kw={'width_ratios': [2.3, 1]})
    fig.subplots_adjust(left=.255, right=.975, bottom=.13, top=.865, wspace=.19)
    y = np.arange(len(rows))
    score_ax.set_yticks(y, [f'{name} (n = {len(frame):,})' for name, frame in rows])
    score_ax.set_ylim(len(rows) - .45, -.65)
    score_ax.set_xlim(xlim)
    grid = np.linspace(*xlim, 1024)
    plotted, support = 0, 0
    for index, (name, frame) in enumerate(rows):
        part = frame.loc[np.isfinite(frame.score)].copy()
        if part.empty:
            score_ax.text(.02, index, 'Not scored', transform=score_ax.get_yaxis_transform(), va='center', fontsize=8.5)
            continue
        values = part.score.to_numpy(float)
        sd = values.std(ddof=1)
        density = gaussian_kde(values, bw_method=.30 / sd)(grid) if sd > 0 else norm.pdf(grid, loc=values[0], scale=.30)
        width = .33 * density / density.max()
        color = '#7DA88C' if index == 0 else ('#A3A8AB' if index < 3 else SAME)
        score_ax.fill_between(grid, index - width, index + width, color=color, alpha=.13, linewidth=0)
        score_ax.plot(grid, index - width, color=color, lw=.55)
        score_ax.plot(grid, index + width, color=color, lw=.55)
        jitter = np.random.default_rng(903 + index).uniform(-.105, .105, len(part))
        scatter = score_ax.scatter(values, index + jitter, s=.65 if index < 3 else 4,
            color=color, alpha=.055 if index < 3 else .55, edgecolors='none', zorder=3)
        assert len(scatter.get_offsets()) == len(part)
        red = part.indel_support.to_numpy(bool)
        score_ax.scatter(values[red], index + jitter[red], s=17, color=IND_WITH_NTC, edgecolors='none', zorder=5)
        plotted += len(part)
        support += int(red.sum())
    score_ax.axvline(0, color='#8F8F8F', lw=.65, zorder=0)
    score_ax.axhline(2.5, color='#D2D2D2', lw=.5)
    score_ax.set_xlabel('Mixscale score (NTC-standardized)')
    score_ax.xaxis.set_major_locator(MaxNLocator(5))
    count_ax.set_xlim(0, bar_limit)
    for index, (name, frame) in enumerate(rows):
        if index < 3:
            count_ax.text(.02, index, 'Reference' if index == 0 else 'Additional control',
                          transform=count_ax.get_yaxis_transform(), va='center', color='#777777', fontsize=8)
            continue
        item = counts.loc[name]
        if item.assessment != 'evaluated':
            count_ax.text(.02, index, 'Not assessed', transform=count_ax.get_yaxis_transform(), va='center', fontsize=8)
            continue
        number = int(item.DE_RNAs)
        count_ax.barh(index, number, height=.52, color=SAME, alpha=.85)
        count_ax.text(number + .025 * bar_limit, index, str(number), ha='left', va='center', fontsize=8.5)
    count_ax.set_xlabel('DE RNAs\n(SCEPTRE, FDR < 0.1)')
    count_ax.xaxis.set_major_locator(MaxNLocator(4, integer=True))
    count_ax.tick_params(axis='y', left=False, labelleft=False)
    for ax in (score_ax, count_ax):
        ax.spines[['top', 'right', 'left']].set_visible(False)
        ax.tick_params(axis='y', length=0)
    fig.suptitle(f'{target} · All chips', fontsize=12, y=.975)
    fig.legend(handles=[Line2D([], [], linestyle='none', marker='o', color=IND_WITH_NTC,
        markeredgewidth=0, markersize=4.5, label='RNA indel support')],
        loc='upper left', bbox_to_anchor=(.255, .955), fontsize=9, frameon=False)
    fig.savefig(destination, dpi=300)
    positions = [list(ax.get_position().bounds) for ax in (score_ax, count_ax)]
    plt.close(fig)
    return {'target': target, 'finite_score_records': plotted, 'plotted_indel_records': support,
            'rows': [{'label': label, 'records': len(frame), 'scored': int(frame.score.notna().sum())} for label, frame in rows],
            'axes': positions, 'sha256': hashlib.sha256(destination.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--labels', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sceptre-stage', default='sceptre')
    args = parser.parse_args()
    cases = load_inputs(args.input, args.labels, args.sceptre_stage)
    summaries = pd.read_csv(args.input / 'mixscale/mixscale_summary.csv').set_index('target_gene_name')
    finite = [frame.score.dropna().to_numpy() for rows, _ in cases.values() for _, frame in rows if frame.score.notna().any()]
    if finite:
        values = np.concatenate(finite)
        margin = max(.3, .025 * np.ptp(values))
        xlim = (float(values.min() - margin), float(values.max() + margin))
    else:
        xlim = (-3, 3)
    maximum = max([int(counts.DE_RNAs.max()) for _, counts in cases.values()] + [1])
    bar_limit = maximum * 1.23
    font = Path('C:/Windows/Fonts/arial.ttf')
    configure(font if font.is_file() else None)
    args.output.mkdir(parents=True, exist_ok=False)
    receipts = [draw(target, rows, counts, xlim, bar_limit, args.output / f'{target}_RNA_response.png',
                     int(summaries.loc[target, 'initial_DE_features']))
                for target, (rows, counts) in cases.items()]
    assert receipts[0]['axes'] == receipts[1]['axes']
    (args.output / 'render_receipt.json').write_text(json.dumps({
        'score_limits': xlim, 'common_density_bandwidth': .30, 'cases': receipts,
        'labels_sha256': hashlib.sha256(args.labels.read_bytes()).hexdigest(),
    }, indent=2) + '\n')


if __name__ == '__main__':
    main()
