"""Plot the observed SCEPTRE NTC calibration without inventing independence."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from render_share_ready_figures import configure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    source = args.input / 'NTC_calibration_tests.csv.gz'
    frame = pd.read_csv(source)
    assert len(frame) == 5000 and frame.grna_target.eq('non-targeting').all()
    observed = np.sort(frame.p_value.to_numpy(float))
    assert np.isfinite(observed).all() and (observed > 0).all() and (observed <= 1).all()
    expected = (np.arange(len(observed)) + .5) / len(observed)
    configure(Path('C:/Windows/Fonts/arial.ttf'))
    fig, ax = plt.subplots(figsize=(4.7, 4.3))
    fig.subplots_adjust(left=.16, right=.95, bottom=.15, top=.80)
    x, y = -np.log10(expected), -np.log10(observed)
    limit = max(x.max(), y.max()) + .1
    ax.plot([0, limit], [0, limit], color='#AAAAAA', lw=.8)
    ax.plot(x[::-1], y[::-1], color='#356A8B', lw=1)
    ax.set(xlim=(0, limit), ylim=(0, limit), xlabel='Expected −log10(P)', ylabel='Observed −log10(P)')
    ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('SCEPTRE · NTC calibration', fontsize=12, y=.97)
    fraction = float(np.mean(observed < .05))
    bh_count = int(frame.significant.fillna(False).sum())
    fig.text(.16, .875, f'5,000 tests · P < 0.05: {fraction:.2%} · BH discoveries: {bh_count}', fontsize=8.5)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300)
    plt.close(fig)
    receipt = {'tests': len(frame), 'fraction_P_below_05': fraction,
               'BH_discoveries_at_10': bh_count,
               'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
               'figure_sha256': hashlib.sha256(args.output.read_bytes()).hexdigest()}
    args.output.with_suffix('.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    main()
