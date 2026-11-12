"""Show the unchanged Mixscale feature gate across declared sensitivity checks."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

from render_share_ready_figures import configure


SCENARIOS = ['NTC_all_guides', 'NTC_guides_ge100', 'intergenic_all_guides',
             'intergenic_guides_ge100', 'no_guide_all_guides']
LABELS = ['NTC · all sgRNAs', 'NTC · sgRNAs with ≥100 records', 'Intergenic · all sgRNAs',
          'Intergenic · sgRNAs with ≥100 records', 'No sgRNA detected · all sgRNAs']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = args.input / 'feature_eligibility.csv'
    frame = pd.read_csv(source)
    assert len(frame) == 10 and not frame.duplicated(['target_gene_name', 'scenario']).any()
    assert frame.required_features.eq(5).all()
    configure(Path('C:/Windows/Fonts/arial.ttf'))
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=.405, right=.96, bottom=.18, top=.72, wspace=.22)
    maximum = max(5, int(frame.DE_features.max())) * 1.25
    for target, ax in zip(['TYMSOS', 'UBA1'], axes):
        part = frame.loc[frame.target_gene_name.eq(target)].set_index('scenario').loc[SCENARIOS]
        y = np.arange(5)
        ax.barh(y, part.DE_features, height=.48, color='#356A8B')
        for index, value in enumerate(part.DE_features):
            ax.text(value + .025 * maximum, index, str(value), va='center', fontsize=10)
        ax.axvline(5, color='#888888', linestyle='--', lw=.8)
        ax.set_xlim(0, maximum)
        ax.set_ylim(4.6, -.7)
        ax.set_title(target, fontsize=11)
        ax.set_xlabel('Eligible DE RNAs')
        ax.set_yticks(y, LABELS)
        ax.xaxis.set_major_locator(MaxNLocator(4, integer=True))
        ax.tick_params(axis='y', length=0, pad=10)
        ax.spines[['top', 'right', 'left']].set_visible(False)
    fig.suptitle('Mixscale · feature eligibility', fontsize=12, y=.96)
    fig.text(.405, .865, 'Dashed line: 5 RNAs required for a score model', fontsize=9)
    fig.text(.405, .07, 'The primary NTC comparison retains all five sgRNAs.', fontsize=9, color='#555555')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300)
    plt.close(fig)
    receipt = {'rows': frame.to_dict('records'), 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
               'figure_sha256': hashlib.sha256(args.output.read_bytes()).hexdigest()}
    args.output.with_suffix('.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
