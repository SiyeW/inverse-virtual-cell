"""Redraw six meeting figures from frozen result tables; no scientific reruns."""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import subprocess
import sys
import textwrap
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors, font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Rectangle
from matplotlib.text import Text
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, norm

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PROJECT = REPO.parent
sys.path.insert(0, str(REPO / 'scripts'))
from plot_mixscale_lower_minimum import load_case

SIZE = (40 / 3, 7.5)
BLUE = '#3F7F9D'
GRAY = '#B8C0C6'
BROWN = '#7A6352'
GREEN = '#AEDBB5'
RED = '#D73936'
DOWN = '#B47C4B'
INK = '#262626'
SOURCES = {
    'old_summary': PROJECT / 'operations/20261026_method_comparison/company_results/scanpy_wilcoxon_pilot_summary.csv',
    'sceptre_counts': PROJECT / 'operations/20261112_two_target_rna/results/sceptre_recovery_r3/guide_DE_counts.csv',
    'sceptre_full': PROJECT / 'operations/20261112_two_target_rna/results/sceptre_recovery_r3/guide_RNA_tests.csv.gz',
    'sceptre_significant': PROJECT / 'operations/20261112_two_target_rna/significant_RNA_details.csv',
    'mixscale_original': PROJECT / 'operations/20261112_two_target_rna/results/mixscale/mixscale_summary.csv',
    'minimum_summary': PROJECT / 'operations/20261113_mixscale_lower_minimum/results/threshold_summary.csv',
    'score_TYMSOS': PROJECT / 'operations/20261113_mixscale_lower_minimum/results/TYMSOS_minimum1_scores.csv.gz',
    'score_UBA1': PROJECT / 'operations/20261113_mixscale_lower_minimum/results/UBA1_minimum1_scores.csv.gz',
    'score_summary': PROJECT / 'reports/20261029_two_figure_briefing/mixscale_lower_minimum_20261113/group_score_summary.csv',
    'score_render': PROJECT / 'reports/20261029_two_figure_briefing/mixscale_lower_minimum_20261113/render_receipt.json',
    'indel_labels': PROJECT / 'operations/20261105_frameshift_RNA_states/coding_rna_record_labels.json',
    'reused_score_loader': REPO / 'scripts/plot_mixscale_lower_minimum.py',
}


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def style():
    font_manager.fontManager.addfont('C:/Windows/Fonts/arial.ttf')
    plt.rcParams.update({
        'font.family': 'Arial', 'font.size': 12, 'axes.labelsize': 13,
        'axes.titlesize': 15, 'axes.titleweight': 'normal', 'axes.linewidth': .8,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.unicode_minus': False,
        'xtick.labelsize': 12, 'ytick.labelsize': 12, 'xtick.major.size': 4,
        'ytick.major.size': 4, 'legend.fontsize': 11.5, 'legend.frameon': False,
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
        'savefig.facecolor': 'white', 'svg.fonttype': 'none', 'pdf.fonttype': 42,
        'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
    })


def sources_manifest():
    return {key: {'path': str(path), 'bytes': path.stat().st_size, 'sha256': sha(path)}
            for key, path in SOURCES.items()}


def load_results():
    data = {'old': pd.read_csv(SOURCES['old_summary']),
            'counts': pd.read_csv(SOURCES['sceptre_counts']),
            'sig': pd.read_csv(SOURCES['sceptre_significant']),
            'original': pd.read_csv(SOURCES['mixscale_original']),
            'minimum': pd.read_csv(SOURCES['minimum_summary']),
            'groups': pd.read_csv(SOURCES['score_summary'])}
    # Read already-computed effects and significance flags, without recalculating P or q.
    full = pd.read_csv(SOURCES['sceptre_full'],
        usecols=['grna_id', 'response_id', 'log_2_fold_change', 'significant'],
        dtype={'significant': 'boolean'}, low_memory=False)
    keep = full.grna_id.str.startswith('TYMSOS-') & full.response_id.isin(
        ['ENSG00000251363', 'ENSG00000289479'])
    repeated = full.loc[keep].copy()
    qmap = data['sig'][['grna_id', 'response_id', 'BH_FDR_check']]
    repeated = repeated.merge(qmap, on=['grna_id', 'response_id'], how='left', validate='one_to_one')
    assert len(repeated) == 10 and int(repeated.significant.sum()) == 4
    data['repeated'] = repeated
    labels = json.loads(SOURCES['indel_labels'].read_bytes())
    root = SOURCES['score_TYMSOS'].parent
    data['cases'] = {target: load_case(root, labels, target) for target in ['TYMSOS', 'UBA1']}
    discrepancy = []
    old_expected = {'TYMSOS': 9, 'UBA1': 11, 'NTC::Non-Targeting-6': 25,
        'NTC::Non-Targeting-40': 24, 'NTC::Non-Targeting-47': 18,
        'NTC::Non-Targeting-56': 12, 'NTC::Non-Targeting-70': 22, 'NTC::Non-Targeting-81': 22}
    old = data['old'].set_index('group')
    for key, expected in old_expected.items():
        if int(old.loc[key, 'nominal_bh05']) != expected:
            discrepancy.append(f'{key}: saved DEG count differs from {expected}')
    expected = {'TYMSOS': ([18, 223, 108, 275, 57], [0, 5, 2, 2, 2]),
                'UBA1': ([16, 86, 99, 123, 214], [1, 0, 1, 0, 2])}
    for target, (sizes, degs) in expected.items():
        counts = data['counts'].set_index('guide_id').loc[[f'{target}-{i}' for i in range(1, 6)]]
        if counts.records.tolist() != sizes or counts.DE_RNAs.tolist() != degs:
            discrepancy.append(f'{target}: SCEPTRE count or record inventory differs')
        rows = data['cases'][target][1]
        stored = data['groups'].loc[data['groups'].target.eq(target)].set_index('group')
        reference = float(stored.loc['NTC', 'NTC_95th_percentile'])
        assert np.isclose(reference, np.quantile(rows[0][1].score, .95))
        for name, frame in rows:
            item = stored.loc[name]
            assert len(frame) == int(item.records)
            assert int(frame.score.gt(reference).sum()) == int(item.above_NTC_95th)
            assert int(frame.indel_support.sum()) == int(item.RNA_indel_support)
            assert int((frame.indel_support & frame.score.gt(reference)).sum()) == int(item.indel_above_NTC_95th)
        assert data['original'].loc[data['original'].target_gene_name.eq(target), 'initial_DE_features'].iloc[0] == 1
        low = data['minimum'].loc[data['minimum'].target_gene_name.eq(target)]
        assert low.available_DE_genes.eq(1).all()
        assert low.loc[low.minimum_DE_genes.gt(1), 'scored_target_records'].eq(0).all()
        assert int(low.loc[low.minimum_DE_genes.eq(1), 'scored_target_records'].iloc[0]) == sum(sizes)
    assert len(data['sig']) == 15 and data['sig'].BH_FDR_check.lt(.1).all()
    assert data['sig'].groupby('grna_id').size().to_dict() == {
        row.guide_id: row.DE_RNAs for row in data['counts'].itertuples() if row.DE_RNAs}
    data['discrepancies'] = discrepancy
    return data


def save(fig, stem, details):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    width, height = fig.canvas.get_width_height()
    outside = []
    for text in fig.findobj(match=Text):
        if not text.get_visible() or not text.get_text().strip():
            continue
        box = text.get_window_extent(renderer)
        if box.x0 < -2 or box.y0 < -2 or box.x1 > width + 2 or box.y1 > height + 2:
            outside.append(text.get_text())
    fig.savefig(HERE / f'{stem}.png', dpi=300)
    fig.savefig(HERE / f'{stem}.svg')
    vector = (HERE / f'{stem}.svg').read_text(encoding='utf-8')
    assert '<text' in vector and '<image' not in vector
    plt.close(fig)
    return {'stem': stem, 'details': details, 'out_of_canvas_text': outside,
            'png_sha256': sha(HERE / f'{stem}.png'), 'svg_sha256': sha(HERE / f'{stem}.svg')}


def fig1(data):
    order = ['TYMSOS', 'UBA1'] + [f'NTC::Non-Targeting-{i}' for i in [6, 40, 47, 56, 70, 81]]
    table = data['old'].set_index('group').loc[order]
    assert table.reference_cells.eq(16952).all()
    fig, ax = plt.subplots(figsize=SIZE)
    fig.subplots_adjust(left=.31, right=.91, top=.79, bottom=.18)
    y = np.arange(8)
    ax.barh(y, table.nominal_bh05, height=.55, color=[BLUE] * 2 + [GRAY] * 6)
    labels = [name.replace('NTC::', 'NTC · ') + f' (n = {int(n):,})' for name, n in zip(order, table.cells)]
    ax.set_yticks(y, labels)
    ax.set_ylim(7.6, -.6)
    ax.set_xlim(0, 30)
    ax.set_xticks(range(0, 31, 5))
    for pos, value in enumerate(table.nominal_bh05):
        ax.text(value + .4, pos, str(int(value)), va='center', fontsize=14)
    ax.axhline(1.5, color='#DDDDDD', lw=.7)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.tick_params(axis='y', length=0, pad=10)
    ax.set_xlabel('Differentially expressed gene/ORF features')
    fig.suptitle('1 · Supplement: conventional DEG pilot and NTC groups', fontsize=19, y=.965)
    fig.text(.31, .875, 'Scanpy Wilcoxon · 1,719 HVGs · BH-adjusted P < 0.05', fontsize=13)
    fig.legend(handles=[Line2D([], [], color=BLUE, lw=8, label='Target groups'),
        Line2D([], [], color=GRAY, lw=8, label='NTC sgRNA groups')],
        loc='upper left', bbox_to_anchor=(.305, .84), ncol=2, fontsize=11.5,
        handlelength=1.2, columnspacing=2)
    fig.text(.31, .065, 'Reference: 16,952 NTC records; the six displayed NTC guides were excluded from this reference.', fontsize=11)
    return save(fig, 'Fig1_conventional_DEG_NTC', {'groups': order, 'DEG_counts': table.nominal_bh05.tolist(),
                'records': table.cells.tolist(), 'reference_records': 16952, 'tested_HVGs': 1719})


def qtext(value):
    return f'{value:.3g}' if value >= .001 else f'{value:.2e}'


def fig2(data):
    order = [f'{t}-{i}' for t in ['TYMSOS', 'UBA1'] for i in range(1, 6)]
    table = data['counts'].set_index('guide_id').loc[order]
    fig = plt.figure(figsize=SIZE)
    ax = fig.add_axes([.19, .18, .30, .62])
    heat = fig.add_axes([.645, .245, .26, .49])
    y = np.arange(10)
    ax.barh(y, table.DE_RNAs, height=.52, color=[BLUE] * 5 + ['#8B9CA9'] * 5)
    ax.set_ylim(9.6, -.6)
    ax.set_yticks(y, [f'{name} (n = {int(n)})' for name, n in zip(order, table.records)])
    ax.set_xlim(0, 6)
    ax.set_xticks(range(7))
    ax.set_xlabel('Significant DEG count')
    ax.tick_params(axis='y', length=0, pad=10)
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.axhline(4.5, color='#DDDDDD', lw=.7)
    for index, value in enumerate(table.DE_RNAs):
        ax.text(value + .12, index, str(int(value)), va='center', fontsize=13)
    genes = ['ENSG00000251363', 'ENSG00000289479']
    rows = data['repeated'].set_index(['grna_id', 'response_id'])
    matrix = np.array([[rows.loc[(f'TYMSOS-{i}', gene), 'log_2_fold_change'] for gene in genes] for i in range(1, 6)])
    limit = max(1, np.ceil(np.abs(matrix).max() * 2) / 2)
    cmap = colors.LinearSegmentedColormap.from_list('direction', [DOWN, '#FFFFFF', BLUE])
    image = heat.pcolormesh(np.arange(3) - .5, np.arange(6) - .5, matrix,
                           vmin=-limit, vmax=limit, cmap=cmap, shading='flat')
    heat.set_xlim(-.5, 1.5);heat.set_ylim(4.5, -.5)
    heat.set_xticks([0, 1], ['LINC02315', 'ENSG00000289479'])
    heat.xaxis.tick_top()
    heat.tick_params(axis='both', length=0, pad=9, labelsize=11.5)
    heat.set_yticks(range(5), [f'sgRNA {i}' for i in range(1, 6)])
    heat.spines[:].set_visible(False)
    for i in range(5):
        for j, gene in enumerate(genes):
            item = rows.loc[(f'TYMSOS-{i + 1}', gene)]
            significant = bool(item.significant)
            value = float(item.log_2_fold_change)
            label = f'{value:+.2f}\n' + (f'q = {qtext(item.BH_FDR_check)}' if significant else 'n.s.')
            text_color = 'white' if abs(value) > .65 * limit else INK
            heat.text(j, i, label, ha='center', va='center', fontsize=11.5, color=text_color)
            heat.add_patch(Rectangle((j - .5, i - .5), 1, 1, fill=False,
                edgecolor='#333333' if significant else 'white', lw=1.3 if significant else .7))
    cax = fig.add_axes([.93, .245, .014, .49])
    colorbar = fig.colorbar(image, cax=cax)
    colorbar.solids.set_rasterized(False)
    colorbar.ax.set_title('log2FC', fontsize=11.5, pad=10)
    colorbar.ax.tick_params(labelsize=11)
    fig.suptitle('2 · SCEPTRE differential expression', fontsize=19, y=.965)
    fig.text(.19, .90, 'Per-sgRNA comparisons with NTC · BH FDR < 0.1', fontsize=13)
    fig.text(.19, .845, 'All five sgRNAs retained per target', fontsize=12.5)
    fig.text(.645, .845, 'TYMSOS · repeated significant genes', fontsize=12.5)
    fig.text(.19, .095, 'NTC: 19,342 records. BH correction covers 273,040 valid guide–feature tests.', fontsize=11)
    fig.text(.19, .052, 'Heatmap: saved log2FC; outlined cells show BH q < 0.1. n.s.: not significant at this threshold.', fontsize=11)
    return save(fig, 'Fig2_SCEPTRE_counts_repeated_genes', {'guide_order': order, 'DEG_counts': table.DE_RNAs.tolist(),
        'heatmap_genes': genes, 'heatmap_log2FC': matrix.tolist(), 'significant_heatmap_cells': 4})


def box(ax, bounds, text, *, fill='#F6F8F9', fontsize=13):
    x, y, w, h = bounds
    ax.add_patch(Rectangle((x, y), w, h, facecolor=fill, edgecolor='#9BA6AD', lw=.85))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fontsize, linespacing=1.35)


def arrow(ax, start, end):
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle='-|>', mutation_scale=13, color='#596A75', lw=1))


def fig3(data):
    fig, ax = plt.subplots(figsize=SIZE)
    ax.set(xlim=(0, 1), ylim=(0, 1));ax.axis('off')
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig.suptitle('3 · Independent DEG-testing and scoring routes', fontsize=19, y=.97)
    box(ax, (.12, .825, .76, .095), 'Saved RNA counts and frozen sgRNA assignments\nTYMSOS: 681 records · UBA1: 538 records · NTC: 19,342 records', fontsize=13)
    arrow(ax, (.35, .825), (.27, .775));arrow(ax, (.65, .825), (.74, .775))
    ax.text(.27, .75, 'SCEPTRE · per sgRNA', ha='center', fontsize=15, color=BLUE)
    ax.text(.74, .75, 'Mixscale · per target', ha='center', fontsize=15, color=BLUE)
    ax.plot([.5, .5], [.13, .76], color='#DDDDDD', lw=.8, ls='--')
    box(ax, (.075, .585, .39, .11), 'One sgRNA group vs NTC\nRaw RNA counts')
    arrow(ax, (.27, .585), (.27, .535))
    box(ax, (.075, .425, .39, .11), 'SCEPTRE testing\nPermutations · BH FDR < 0.1')
    arrow(ax, (.27, .425), (.27, .365))
    box(ax, (.075, .235, .39, .13), 'Output: per-sgRNA DEG lists\nTYMSOS counts: 0, 5, 2, 2, 2\nUBA1 counts: 1, 0, 1, 0, 2', fontsize=12.5)
    ax.text(.27, .16, 'Computed independently of the\nMixscale DEG selection', ha='center', va='center', fontsize=12)
    box(ax, (.545, .59, .39, .105), 'Pool all five sgRNAs for each target\nSeurat DEG screen vs NTC')
    arrow(ax, (.74, .59), (.74, .555))
    box(ax, (.545, .425, .39, .13), 'Wilcoxon · Bonferroni-adjusted P < 0.05\nTYMSOS: LINC02315 (1 gene)\nUBA1: ZFP42 (1 gene)', fontsize=12)
    arrow(ax, (.74, .425), (.74, .39))
    box(ax, (.545, .28, .39, .11), 'Original min.de.genes = 5: not scored\nLowered to 1: overall score available', fontsize=12.5)
    arrow(ax, (.74, .28), (.74, .245))
    box(ax, (.545, .13, .39, .115), 'Exploratory single-gene score\n10 shared PCs + 20-nearest-NTC background\nPCA fit: 81,201 input records', fontsize=11.5)
    ax.text(.06, .055, 'The minimum DEG count controls entry to scoring; it is not a measure of overall experiment quality.', fontsize=11.5)
    return save(fig, 'Fig3_Mixscale_scoring_limits', {'independent_routes': True,
        'available_genes_per_target': 1, 'original_minimum': 5, 'exploratory_minimum': 1})


def score_figure(data, target, number):
    marker, rows = data['cases'][target]
    stored = data['groups'].loc[data['groups'].target.eq(target)].set_index('group')
    original = json.loads(SOURCES['score_render'].read_text())
    xlim = original['score_limits'];tail_limit = original['percent_limits'][1]
    fig, (ax, tail) = plt.subplots(1, 2, figsize=SIZE, sharey=True,
        gridspec_kw={'width_ratios': [2.65, 1]})
    fig.subplots_adjust(left=.26, right=.965, top=.79, bottom=.17, wspace=.25)
    ax.set_yticks(range(8), [f'{name} (n = {len(frame):,})' for name, frame in rows])
    ax.set_ylim(7.55, -.6);ax.set_xlim(xlim);tail.set_xlim(0, tail_limit)
    reference = float(stored.loc['NTC', 'NTC_95th_percentile'])
    ax.axvline(0, color='#DDDDDD', lw=.6, zorder=0)
    ax.axvline(reference, color='#777777', lw=.9, ls='--', zorder=1)
    grid = np.linspace(*xlim, 900)
    total, red_total = 0, 0
    for i, (name, frame) in enumerate(rows):
        values = frame.score.to_numpy(float);sd = values.std(ddof=1)
        density = gaussian_kde(values, bw_method=.25 / sd)(grid) if sd > 0 else norm.pdf(grid, values[0], .25)
        width = .3 * density / density.max();support = np.flatnonzero(density > .002 * density.max())
        extent = slice(support[0], support[-1] + 1)
        color = GREEN if i == 0 else (GRAY if i < 3 else BROWN)
        fill = '#E4F0E6' if i == 0 else ('#F0F1F2' if i < 3 else '#EEE8E3')
        ax.fill_between(grid[extent], i - width[extent], i + width[extent], color=fill, linewidth=0)
        ax.plot(grid[extent], i - width[extent], color=color, lw=.7)
        ax.plot(grid[extent], i + width[extent], color=color, lw=.7)
        jitter = np.random.default_rng(1130 + i).uniform(-.085, .085, len(frame))
        if i < 3:
            points = ax.scatter(values, i + jitter, s=.85, color=color, edgecolors='none', zorder=2)
        else:
            points = ax.scatter(values, i + jitter, s=10, color=color, marker='x', linewidths=.65, zorder=3)
        assert len(points.get_offsets()) == len(frame)
        red = frame.indel_support.to_numpy(bool)
        highlights = ax.scatter(values[red], i + jitter[red], s=25, color=RED, edgecolors='none', zorder=5)
        assert len(highlights.get_offsets()) == int(red.sum())
        total += len(frame);red_total += int(red.sum())
        record = stored.loc[name];pct = float(record.percent_above_NTC_95th)
        tail.barh(i, pct, height=.43, color=color)
        tail.text(pct + .035 * tail_limit, i, f'{int(record.above_NTC_95th)}/{len(frame):,}\n({pct:.1f}%)',
            va='center', fontsize=11.5)
    for item in [ax, tail]:
        item.spines[['top', 'right', 'left']].set_visible(False)
        item.tick_params(axis='y', length=0)
        item.axhline(2.5, color='#D8D8D8', lw=.7)
    tail.tick_params(axis='y', labelleft=False, left=False)
    ax.xaxis.set_major_locator(MaxNLocator(6));tail.set_xticks(range(0, 41, 10))
    ax.set_xlabel('Single-gene score\n(NTC mean = 0, SD = 1)')
    tail.set_xlabel('Above NTC 95th percentile\n(% of group records)')
    fig.suptitle(f'{number} · {target} · Exploratory single-gene score', fontsize=19, y=.97)
    fig.text(.26, .91, f'Scoring gene: {marker} · min.de.genes = 1 (original setting: 5)', fontsize=13)
    handles = [Line2D([], [], color=BROWN, marker='x', ls='none', markersize=5, label='Targeting records'),
        Line2D([], [], color=RED, marker='o', ls='none', markersize=5, label='RNA indel support'),
        Line2D([], [], color='#777777', ls='--', lw=.9, label='NTC 95th percentile')]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.252, .875), ncol=3,
        fontsize=11.5, columnspacing=1.6, handletextpad=.6)
    fig.text(.26, .04, 'Same-dataset exploration. The percentile line is a descriptive reference.', fontsize=11.5, color='#555555')
    result = save(fig, f'Fig{number}_{target}_single_gene_scores', {'target': target, 'marker': marker,
        'all_records_plotted': total, 'indel_support_plotted': red_total, 'NTC_95th_percentile': reference,
        'axes': [list(item.get_position().bounds) for item in [ax, tail]], 'score_limits': xlim})
    return result


def indel_table(data):
    rows = []
    expected = {'TYMSOS': (0, 9, 73, 672), 'UBA1': (3, 31, 64, 507)}
    for target in ['TYMSOS', 'UBA1']:
        part = data['groups'].loc[data['groups'].target.eq(target) & data['groups'].group.str.startswith(target + '-')]
        a = int(part.indel_above_NTC_95th.sum());n = int(part.RNA_indel_support.sum())
        b = int(part.above_NTC_95th.sum()) - a;m = int(part.records.sum()) - n
        assert (a, n, b, m) == expected[target]
        rows += [{'target': target, 'category': 'RNA indel support', 'above_reference': a, 'records': n, 'percent': a / n * 100},
                 {'target': target, 'category': 'No observed RNA indel support', 'above_reference': b, 'records': m, 'percent': b / m * 100}]
    return pd.DataFrame(rows)


def fig6(data):
    table = indel_table(data)
    fig, ax = plt.subplots(figsize=SIZE)
    fig.subplots_adjust(left=.19, right=.84, top=.79, bottom=.22)
    positions = [0, .62, 1.9, 2.52]
    for i, row in enumerate(table.itertuples()):
        color = RED if i % 2 == 0 else GRAY
        ax.barh(positions[i], row.percent, height=.43, color=color)
        ax.text(row.percent + .3, positions[i], f'{row.percent:.2f}%  ({row.above_reference}/{row.records})',
                fontsize=14, va='center')
        if row.percent == 0:
            ax.plot([0], [positions[i]], 'o', ms=4, mfc='white', mec=RED, clip_on=False)
    ax.set_yticks([.31, 2.21], ['TYMSOS', 'UBA1'])
    ax.tick_params(axis='y', length=0, pad=18, labelsize=15)
    ax.set_ylim(3.1, -.65);ax.set_xlim(0, 16)
    ax.set_xticks(range(0, 17, 4))
    ax.set_xlabel('Records above the target-specific NTC 95th percentile (%)')
    ax.spines[['top', 'right', 'left']].set_visible(False)
    ax.axhline(1.25, color='#DDDDDD', lw=.7)
    handles = [Line2D([], [], color=RED, lw=9, label='RNA indel support'),
        Line2D([], [], color=GRAY, lw=9, label='No observed RNA indel support')]
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.18, .845), ncol=2, fontsize=12.5,
               handlelength=1.2, columnspacing=2)
    fig.suptitle('6 · RNA indel support and score upper tails', fontsize=19, y=.965)
    fig.text(.19, .90, 'TYMSOS: LINC02315 score · UBA1: ZFP42 score · min.de.genes = 1', fontsize=13)
    fig.text(.19, .115, 'RNA-indel-supported groups are small: n = 9 and n = 31.', fontsize=12)
    fig.text(.19, .065, 'No observed support includes records without locus coverage; it does not establish an unedited genotype.', fontsize=11.5)
    return save(fig, 'Fig6_indel_support_score_upper_tail', {'rows': table.to_dict('records')}), table


def supplement(data):
    fig, ax = plt.subplots(figsize=SIZE)
    ax.axis('off');fig.subplots_adjust(left=.04, right=.96, top=.86, bottom=.14)
    rows = []
    for row in data['sig'].itertuples():
        rows.append([row.grna_id, row.gene_name, f'{row.log_2_fold_change:+.3f}', f'{row.BH_FDR_check:.3g}'])
    table = ax.table(cellText=rows, colLabels=['sgRNA', 'Gene / ORF feature', 'log2FC', 'BH-adjusted P'],
        colWidths=[.16, .54, .12, .18], cellLoc='left', colLoc='left', bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False);table.set_fontsize(11.3)
    for (i, j), cell in table.get_celld().items():
        cell.PAD = .10
        cell.set_facecolor('white');cell.set_edgecolor('#E0E0E0');cell.set_linewidth(.45)
        if i == 0:
            cell.set_text_props(weight='bold');cell.set_edgecolor('#777777');cell.set_linewidth(.8)
        elif j == 2:
            cell.get_text().set_color(BLUE if data['sig'].iloc[i - 1].log_2_fold_change > 0 else DOWN)
    fig.suptitle('Supplement · All 15 significant SCEPTRE associations', fontsize=19, y=.965)
    fig.text(.04, .905, 'BH FDR < 0.1 across 273,040 valid guide–feature tests; all feature identifiers are retained.', fontsize=12)
    fig.text(.04, .065, 'Positive log2FC: increased expression; negative: decreased expression. Full raw P and detection counts are in the CSV.', fontsize=11.3)
    return save(fig, 'S1_SCEPTRE_all15', {'rows': len(rows)})


def documents(data, checks, source_info):
    captions = {
      'Fig1_conventional_DEG_NTC': ('图1（补充图）｜旧普通DEG试跑及NTC对照',
        '读取10月26日Scanpy Wilcoxon试跑的保存汇总，显示TYMSOS、UBA1及六条单独比较的NTC sgRNA。检验范围为1,719个HVG，按每个比较组的BH校正后P<0.05统计，未另加logFC过滤；共同参照为其余16,952条NTC记录。TYMSOS、UBA1分别检出9、11项，六个NTC组为12～25项。对照组间的检出不能直接认定为假阳性；本图为同批合并条码记录的描述。'),
      'Fig2_SCEPTRE_counts_repeated_genes': ('图2｜逐sgRNA的SCEPTRE差异表达结果',
        '左侧显示两个靶点全部五条sgRNA的显著基因/ORF数量及记录数，以19,342条NTC为参照，在273,040个有效guide—表达条目检验中统一采用BH FDR<0.1。右侧读取保存的log₂FC和显著性标记：LINC02315在TYMSOS-2、4中上调，ENSG00000289479在TYMSOS-3、5中下调；描边格显示已有BH校正P，n.s.表示未达该阈值。少量跨guide重复变化不能直接证明DNA编辑或功能后果。'),
      'Fig3_Mixscale_scoring_limits': ('图3｜两条独立计算路线与Mixscale评分条件',
        'SCEPTRE按单条sgRNA检验；Mixscale另将同靶点五条sgRNA合并，使用Seurat Wilcoxon及Bonferroni校正后P<0.05重新选择评分基因，两条路线没有把DEG清单直接相互传递。TYMSOS仅得到LINC02315、UBA1仅得到ZFP42，均不足原最低5个；将min.de.genes降到1后整体投影可计算，实际为单基因探索评分。最低基因数是评分准入设置，不是实验质量标准；单基因的留一辅助计算不可用。'),
      'Fig4_TYMSOS_single_gene_scores': ('图4｜TYMSOS的单基因探索评分',
        '复用11月13日全部逐记录分数，以LINC02315为唯一评分基因，保留681条靶向记录、19,342条NTC、9,094条无sgRNA检出和195条intergenic记录。每个点或叉代表一条记录，红点为9条已有RNA indel支持；外形仅表示平滑后的分布。虚线固定为原NTC实际第95百分位1.736215，右侧列出超过该线的数量/本组总数及比例；TYMSOS-2、4分别为28/223（12.6%）、35/275（12.7%）。参照线只用于描述，同批选基因及估计方向不构成独立分类验证。'),
      'Fig5_UBA1_single_gene_scores': ('图5｜UBA1的单基因探索评分',
        '复用11月13日ZFP42单基因分数，保留538条靶向记录及与图4相同的三类对照，红点为31条已有RNA indel支持。虚线沿用NTC实际第95百分位1.633113，右侧列出超过数量和比例。UBA1-5为58/214（27.1%），但intergenic对照也为39/195（20.0%）；高分并非UBA1靶向记录独有。单基因、同批数据内的探索不能直接解释为已验证的特异扰动响应。'),
      'Fig6_indel_support_score_upper_tail': ('图6｜RNA indel支持与评分上尾的对应',
        '将图4、5中超过各自NTC第95百分位的记录，与现有RNA indel支持标签交叉统计。TYMSOS有支持者0/9（0%），未观察到支持者73/672（10.86%）；UBA1分别为3/31（9.68%）与64/507（12.62%）。图中同时给出绝对数量和比例，不新增显著性检验。支持组样本量仅9和31；未观察到支持包含缺乏位点覆盖者，不能标作未编辑细胞，也不能据此估计DNA编辑率。'),
      'S1_SCEPTRE_all15': ('补充表S1｜全部15项SCEPTRE显著关联',
        '完整列出已有15项显著guide—基因/ORF关联，log₂FC沿用软件输出，BH校正P读取既有复核表，不重新计算。正负值表示相对NTC的升降；所有原始特征标识保留。随附CSV还保留原P值和两组检出记录数。该表不将显著表达关联等同于功能KO。'),
    }
    mapping = {
      'Fig1_conventional_DEG_NTC': ('fig1', ['old_summary']),
      'Fig2_SCEPTRE_counts_repeated_genes': ('fig2', ['sceptre_counts', 'sceptre_full', 'sceptre_significant']),
      'Fig3_Mixscale_scoring_limits': ('fig3', ['sceptre_counts', 'mixscale_original', 'minimum_summary']),
      'Fig4_TYMSOS_single_gene_scores': ('score_figure', ['score_TYMSOS', 'score_summary', 'score_render', 'indel_labels', 'reused_score_loader']),
      'Fig5_UBA1_single_gene_scores': ('score_figure', ['score_UBA1', 'score_summary', 'score_render', 'indel_labels', 'reused_score_loader']),
      'Fig6_indel_support_score_upper_tail': ('fig6', ['score_summary']),
      'S1_SCEPTRE_all15': ('supplement', ['sceptre_significant']),
    }
    with (HERE / '绘图索引.csv').open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['figure', 'title', 'data_tables', 'script', 'function', 'PNG', 'SVG'])
        writer.writeheader()
        for stem, (title, _) in captions.items():
            function, keys = mapping[stem]
            writer.writerow({'figure': stem, 'title': title, 'data_tables': ' | '.join(str(SOURCES[key]) for key in keys),
                'script': str(Path(__file__).resolve()), 'function': function, 'PNG': str(HERE / f'{stem}.png'), 'SVG': str(HERE / f'{stem}.svg')})
    (HERE / '中文图注.txt').write_text('\n\n'.join(title + '\n' + caption for title, caption in captions.values()) + '\n', encoding='utf-8')
    sections = []
    for stem, (title, caption) in captions.items():
        function, keys = mapping[stem]
        sources = '<br>'.join(html.escape(str(SOURCES[key])) for key in keys)
        backup_link = ' · <a href="S1_SCEPTRE_all15.csv">完整结果 CSV</a>' if stem == 'S1_SCEPTRE_all15' else ''
        sections.append(f'<section id="{stem}"><h2>{html.escape(title)}</h2><p><a href="{stem}.png">高分辨率 PNG</a> · <a href="{stem}.svg">可编辑 SVG</a>{backup_link}</p><img src="{stem}.png" alt="{html.escape(title)}"><p>{html.escape(caption)}</p><details><summary>数据与绘图函数</summary><p class="path">{sources}</p><p>make_figures.py · {function}</p></details></section>')
    page = '<!doctype html><html lang="zh-CN"><meta charset="utf-8"><title>Perturb-seq组会绘图索引</title><style>body{max-width:1300px;margin:32px auto;padding:0 24px;font:16px/1.65 Arial,"Microsoft YaHei",sans-serif;color:#262626;background:white}h1{font-size:25px}h2{font-size:20px}section{margin:42px 0 64px}img{width:100%;height:auto;border:1px solid #eee}a{color:#356a8b}.path{font-size:13px;overflow-wrap:anywhere}details{font-size:14px}p{max-width:1150px}</style><h1>2026-11-14 Perturb-seq组会图件</h1><p>仅使用已有输出重绘：六张主图、一张15项结果补充表；PNG为300 dpi，SVG保留可编辑文字和矢量图形。未制作PPT，未重跑科学分析，未覆盖旧结果。</p><p>数据单位均为合并条码记录；真实捕获身份未核实。数字与最近技术交接核对一致。旧口头“NTC 18～25项”的概括不完整，本次按原表12～25项展示。</p>' + ''.join(sections) + '</html>'
    (HERE / 'index.html').write_text(page, encoding='utf-8')
    (HERE / '数值核对.txt').write_text('与最近技术交接及本次用户给定数量核对：未发现数量不一致。\n图1使用nominal_bh05列，不使用附加abs(logFC)>0.25的列。\n旧口头“NTC 18～25项”范围不完整，原保存结果为12～25项；旧文件未改。\n本次没有重新进行DEG检验、P/FDR校正、PCA、评分或方法试验。\n', encoding='utf-8')
    payload = SOURCES['sceptre_significant'].read_bytes()
    backup = HERE / 'S1_SCEPTRE_all15.csv';backup.write_bytes(payload)
    assert sha(backup) == source_info['sceptre_significant']['sha256']
    data['repeated'].to_csv(HERE / 'Fig2_repeated_gene_values.csv', index=False)
    indel_table(data).to_csv(HERE / 'Fig6_cross_counts.csv', index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--only', nargs='*', choices=['1', '2', '3', '4', '5', '6', 'S1'])
    args = parser.parse_args()
    source_info = sources_manifest();data = load_results()
    if data['discrepancies']:
        (HERE / '差异清单.txt').write_text('\n'.join(data['discrepancies']) + '\n', encoding='utf-8')
        raise ValueError('Saved results differ; discrepancy list written. No figures produced.')
    if args.check_only:
        print(json.dumps({'discrepancies': [], 'significant_associations': len(data['sig']),
            'repeated_gene_cells': len(data['repeated']), 'indel_cross_counts': indel_table(data).to_dict('records')}, indent=2))
        return
    style();checks=[]
    selected = args.only or ['1', '2', '3', '4', '5', '6', 'S1']
    for name in selected:
        if name == '1': result = fig1(data)
        elif name == '2': result = fig2(data)
        elif name == '3': result = fig3(data)
        elif name in ['4', '5']: result = score_figure(data, 'TYMSOS' if name == '4' else 'UBA1', name)
        elif name == '6': result, _ = fig6(data)
        else: result = supplement(data)
        checks.append(result);print('RENDERED',result['stem'],'outside',result['out_of_canvas_text'],flush=True)
    existing = HERE / 'render_checks.json'
    prior = json.loads(existing.read_text()) if existing.exists() else []
    names = {check['stem'] for check in checks}
    checks = [check for check in prior if check['stem'] not in names] + checks
    (HERE / 'render_checks.json').write_text(json.dumps(checks, indent=2) + '\n')
    documents(data, checks, source_info)
    after = sources_manifest();assert source_info == after, 'An existing input changed during rendering'
    code_commit = subprocess.check_output(['git', '-C', str(REPO), 'rev-parse', 'HEAD'], text=True).strip()
    reference_decks = json.loads((HERE / 'reference_style.json').read_text())
    assert all(sha(item['source']) == item['sha256'] for item in reference_decks)
    manifest = {'created_at': datetime.now(timezone.utc).isoformat(), 'code_commit': code_commit,
        'script': str(Path(__file__).resolve()), 'script_sha256': sha(__file__), 'sources': source_info,
        'source_files_unchanged': True, 'discrepancies_against_current_handoff': [],
        'reference_decks': reference_decks,
        'figures': checks, 'PPT_created': False, 'scientific_analysis_rerun': False}
    (HERE / 'source_manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
    if all((HERE / (stem + '.png')).exists() for stem in [
        'Fig1_conventional_DEG_NTC', 'Fig2_SCEPTRE_counts_repeated_genes', 'Fig3_Mixscale_scoring_limits',
        'Fig4_TYMSOS_single_gene_scores', 'Fig5_UBA1_single_gene_scores',
        'Fig6_indel_support_score_upper_tail', 'S1_SCEPTRE_all15']):
        from PIL import Image, ImageOps, ImageDraw
        board = Image.new('RGB', (1800, 1560), '#E8E8E8')
        stems = ['Fig1_conventional_DEG_NTC', 'Fig2_SCEPTRE_counts_repeated_genes', 'Fig3_Mixscale_scoring_limits',
                 'Fig4_TYMSOS_single_gene_scores', 'Fig5_UBA1_single_gene_scores', 'Fig6_indel_support_score_upper_tail']
        for i, stem in enumerate(stems):
            with Image.open(HERE / (stem + '.png')) as source:
                thumb = ImageOps.contain(source.convert('RGB'), (890, 500))
            board.paste(thumb, ((i % 2) * 900 + 5, (i // 2) * 520 + 10))
        # Fixed thumbnail width preserves the same 16:9 ratio for every panel.
        board.save(HERE / 'all_figures_preview.png')


if __name__ == '__main__':
    main()
