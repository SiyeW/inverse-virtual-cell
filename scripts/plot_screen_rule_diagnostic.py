"""Explain a saved guide-selection rule; never rerun or validate the RNA tests."""
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_saved(folder):
    folder = Path(folder)
    target = pd.read_csv(folder / 'guide_screen_summary.csv')
    ntc = pd.read_csv(folder / 'NTC_guide_summary.csv')
    energy = pd.read_csv(folder / 'energy_by_guide.csv')
    summary = pd.read_csv(folder / 'screen_summary.csv').iloc[0]
    assert target.guide.is_unique and ntc.guide.is_unique and energy.guide.is_unique
    assert not set(target.guide) & set(ntc.guide)
    assert len(target) == int(summary.eligible_targeting_guides) == 368
    assert len(ntc) == int(summary.NTC_guides) == 97
    assert int(energy.records.sum()) == int(summary.records) == 12591
    aligned = energy.set_index('guide').loc[target.guide]
    assert np.array_equal(target.records.to_numpy(), aligned.records.to_numpy())
    assert np.array_equal(target.role.to_numpy(), aligned.role.to_numpy())
    fields = ['guide', 'guide_q_BY', 'RNA_association_candidate']
    tests = pd.concat([target[fields], ntc[fields]], ignore_index=True)
    assert tests.guide_q_BY.between(0, 1, inclusive='right').all()
    assert tests.RNA_association_candidate.isin([True, False]).all()
    assert np.array_equal(tests.RNA_association_candidate.to_numpy(), tests.guide_q_BY.lt(.05).to_numpy())
    frame = tests.merge(energy, on='guide', how='left', validate='one_to_one')
    assert len(frame) == 465 and frame.notna().all().all()
    assert frame.records.ge(10).all()
    assert set(frame.role) == {'targeting', 'positive_control', 'non_targeting_control'}
    assert set(frame.loc[frame.role.eq('non_targeting_control'), 'guide']) == set(ntc.guide)
    assert frame.loc[~frame.role.eq('non_targeting_control'), 'reference_records'].eq(3769).all()
    nc = frame.loc[frame.role.eq('non_targeting_control')]
    assert nc.records.sum() == int(summary.NTC_records) == 3769
    assert np.array_equal(nc.reference_records, 3769 - nc.records)
    assert int(frame.loc[~frame.role.eq('non_targeting_control'), 'RNA_association_candidate'].sum()) == int(summary.RNA_association_candidates) == 6
    assert int(nc.RNA_association_candidate.sum()) == int(summary.NTC_guides_passing_same_rule) == 5
    untested = energy.loc[~energy.guide.isin(frame.guide)]
    assert len(untested) == 1 and untested.role.eq('intergenic_control').all()
    return frame, untested


def draw(folder, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import ScalarFormatter

    folder, output = Path(folder), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    frame, untested = load_saved(folder)
    groups = [
        ('targeting', '#bd8050', 13, .52, 'Other targeting (n = 351)'),
        ('positive_control', '#e60012', 22, .92, 'Library positive controls (n = 17)'),
        ('non_targeting_control', '#377394', 17, .8, 'NTC (n = 97)'),
    ]
    assert [int(frame.role.eq(x[0]).sum()) for x in groups] == [351, 17, 97]
    selected = frame.RNA_association_candidate
    plt.rcParams.update({'font.sans-serif': ['Arial', 'DejaVu Sans'], 'font.size': 10,
                         'axes.unicode_minus': False, 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.8))
    fig.subplots_adjust(left=.085, right=.975, bottom=.18, top=.80, wspace=.32)
    fig.suptitle('First chip: E-distance and the saved expression-test selection', y=.98, fontsize=12)
    handles = []
    for role, color, size, alpha, label in groups:
        rows = frame.loc[frame.role.eq(role)]
        for ax, value in zip(axes, ['energy_unbiased', 'guide_q_BY']):
            ax.scatter(rows.records, rows[value], s=size, c=color, alpha=alpha,
                       edgecolors='white', linewidths=.2, zorder=2)
        handles.append(Line2D([], [], marker='o', linestyle='none', markersize=4.5,
                              color=color, label=label))
    for ax, value in zip(axes, ['energy_unbiased', 'guide_q_BY']):
        rows = frame.loc[selected]
        ax.scatter(rows.records, rows[value], s=64, facecolors='none', edgecolors='#222222',
                   linewidths=.9, zorder=4)
        ax.set_xscale('log')
        ax.set_xticks([10, 20, 50, 100])
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.minorticks_off()
        ax.set_xlim(9.3, float(frame.records.max()) * 1.12)
        ax.set_xlabel('Cell barcodes per sgRNA', labelpad=7)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(width=.7, length=3.5)
    handles.append(Line2D([], [], marker='o', linestyle='none', markerfacecolor='none',
                          markeredgecolor='#222222', markersize=7,
                          label='Circled: 6 targeting + 5 NTC'))
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(.53, .95), ncol=2,
               frameon=False, fontsize=8.2, columnspacing=2, handletextpad=.5)
    axes[0].axhline(0, color='#999999', linewidth=.6, zorder=1)
    axes[0].set_ylabel('E-distance estimate from NTC', labelpad=6)
    axes[0].set_title('A  Selected sgRNAs in the E-distance view', loc='left', fontsize=10, pad=10)
    axes[1].set_yscale('log')
    axes[1].set_ylim(1e-16, 2)
    axes[1].set_yticks([1e-16, 1e-12, 1e-8, 1e-4, .05, 1])
    axes[1].set_yticklabels([r'$10^{-16}$', r'$10^{-12}$', r'$10^{-8}$', r'$10^{-4}$', '0.05', '1'])
    axes[1].minorticks_off()
    axes[1].axhline(.05, color='#555555', linestyle='--', linewidth=.85, zorder=1)
    axes[1].annotate('P = 0.05', xy=(103, .05), xytext=(0, -5), textcoords='offset points',
                     ha='right', va='top', fontsize=8, color='#444444')
    axes[1].set_ylabel('BY-adjusted P-value (log scale)', labelpad=6)
    axes[1].set_title('B  Actual rule: adjusted P < 0.05', loc='left', fontsize=10, pad=10)
    fig.text(.085, .035, 'All 465 tested sgRNAs are shown. Circles mark all 11 selections; no E-distance cutoff was used.',
             fontsize=8, color='#555555', ha='left')
    basename = '首芯片_E-distance与校正P值筛选规则'
    for ext in ['png', 'svg']:
        fig.savefig(output / f'{basename}.{ext}', dpi=300, facecolor='white')
    plt.close(fig)
    frame.to_csv(output / 'selection_rule_display.csv', index=False)
    proof = dict(
        dataset='first chip, not company merged', tested_guides=len(frame),
        targeting_guides_including_design_positive=368, NTC_guides=97,
        selected_targeting_including_design_positive=6, selected_NTC=5,
        selected_design_positive=1, all_selected_guides_circled=True,
        all_tested_guides_shown=True, threshold_axis='saved BY-adjusted guide P-value', threshold=.05,
        below_threshold_is_selected=True, E_distance_threshold_used=False,
        selected_ids=frame.loc[selected, 'guide'].tolist(),
        intergenic_not_expression_tested=untested[['guide', 'role', 'records']].to_dict('records'),
        all_input_files_unchanged=True, RNA_tests_rerun=False,
        input_SHA256={name: sha(folder/name) for name in
            ['guide_screen_summary.csv', 'NTC_guide_summary.csv', 'energy_by_guide.csv', 'screen_summary.csv']},
        status='visual explanation of the previous exploratory rule, not validation of a reliable selector')
    (output/'plot_proof.json').write_text(json.dumps(proof, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({k: proof[k] for k in ['dataset', 'tested_guides', 'selected_targeting_including_design_positive',
        'selected_NTC', 'all_selected_guides_circled', 'RNA_tests_rerun']}, ensure_ascii=False))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--folder', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    draw(args.folder, args.output)
