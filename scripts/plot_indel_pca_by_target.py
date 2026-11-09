"""Plot each observed insertion or deletion on the accepted PCA for its target."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np

from plot_coding_rna_states import target_display
from render_share_ready_figures import assert_drawn_points, configure


COHORT = {'company': 'Company merged', 'chip01': 'Chip 1'}
OTHER = '#C0C6CB'
SAME = '#626D77'
IND = '#A6403D'


def draw_case(row, arrays, variants, destination):
    xy = arrays['target_PCA'][:, :2]
    variance = arrays['variance_ratio']
    guides = arrays['guide'].astype(str)
    barcodes = arrays['cell_barcode'].astype(str)
    assert arrays['PCA_fit_status'].item() == 'target_fit_NTC_projection'
    assert xy.shape == (int(row['RNA_records']), 2) and np.isfinite(xy).all()
    assert len(set(guides)) == int(row['sgRNAs'])
    lookup = {key: index for index, key in enumerate(zip(guides, barcodes))}
    assert len(lookup) == len(xy)

    panels = variants or [None]
    columns = min(len(panels), 3 if len(panels) <= 9 else 4)
    rows = (len(panels) + columns - 1) // columns
    width, height = 3.25 * columns + .25, 3.1 * rows + 1.1
    fig, grid = plt.subplots(rows, columns, figsize=(width, height),
                             sharex=True, sharey=True, squeeze=False)
    fig.subplots_adjust(left=.58 / width, right=.98, top=1 - .76 / height,
                        bottom=.80 / height, wspace=.21, hspace=.31)
    low, high = xy.min(axis=0), xy.max(axis=0)
    center = (low + high) / 2
    span = .55 * np.ptp(xy, axis=0).max()
    assert span > 0
    annotations = []
    for index, (ax, event) in enumerate(zip(grid.flat, panels)):
        if event is None:
            ax.scatter(*xy.T, s=6, color=OTHER, alpha=.67, edgecolors='none')
            ax.set_title('RNA indels observed = 0', fontsize=9)
        else:
            support = [cell for cell in event['cell_support'] if cell['alt'] > 0]
            assert sum(cell['alt'] for cell in support) == event['guide_alt']
            highlighted = [lookup[(event['guide'], cell['cell_barcode'])] for cell in support]
            assert highlighted and len(set(highlighted)) == len(highlighted)
            same = guides == event['guide']
            same[highlighted] = False
            for group, size, color, alpha in ((xy[guides != event['guide']], 6, OTHER, .62),
                                              (xy[same], 7, SAME, .73),
                                              (xy[highlighted], 19, IND, 1)):
                ax.scatter(*group.T, s=size, color=color, alpha=alpha, edgecolors='none')
            kind = 'deletion' if event['kind'] == 'DEL' else 'insertion'
            number = event['guide'].rsplit('-', 1)[-1]
            location = f"{event['chrom']}:{event['pos0'] + 1:,}"
            length = abs(event['length_change'])
            change = f'{length} bp {kind}'
            if event['kind'] == 'INS' and length <= 4 and event['alt'].startswith(event['ref']):
                change = f"{change} (+{event['alt'][len(event['ref']):]})"
            category = event['category']
            ax.set_title(f'sgRNA {number} · {change} (n = {len(support)})\n{location} · {category}', fontsize=9)
            annotations.append({key: event[key] for key in ('guide', 'chrom', 'pos0', 'ref', 'alt', 'category')}
                               | {'barcode_records': len(support), 'RNA_UMI_families': event['guide_alt']})
        assert_drawn_points(ax, xy)
        ax.set(xlim=(center[0] - span, center[0] + span), ylim=(center[1] - span, center[1] + span),
               aspect='equal')
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        if index + columns >= len(panels):
            ax.tick_params(axis='x', labelbottom=True)
    for ax in list(grid.flat)[len(panels):]:
        ax.set_visible(False)
    fig.suptitle(f"{target_display(row)} · {COHORT[row['dataset']]} ({len(xy)} records)",
                 fontsize=10.5, y=1 - .08 / height)
    fig.supxlabel(f'PC 1 ({100 * variance[0]:.1f}%)', fontsize=10, y=.34 / height)
    fig.supylabel(f'PC 2 ({100 * variance[1]:.1f}%)', fontsize=10, x=.035 / width)
    legend = [(OTHER, 'Target records')]
    if variants:
        legend = [(OTHER, 'Other sgRNAs'), (SAME, 'Same sgRNA'), (IND, 'Supports shown indel')]
    handles = [Line2D([], [], linestyle='none', marker='o', color=color, markersize=4, label=label)
               for color, label in legend]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.5, .02 / height), ncol=len(handles),
               fontsize=8, handletextpad=.4, columnspacing=1.4)
    fig.savefig(destination, dpi=300, bbox_inches='tight', pad_inches=.06)
    plt.close(fig)
    return annotations


def gallery(path, entries):
    fields = ('id', 'name', 'dataset', 'png', 'indels')
    data = json.dumps([{key: r[key] for key in fields} for r in entries], ensure_ascii=False).replace('</', '<\\/')
    path.write_text('''<!doctype html><html lang="en"><meta charset="utf-8">
<title>RNA indels by target</title>
<style>
body{font:15px Arial,sans-serif;color:#242424;margin:24px;background:white}
main{max-width:1320px;margin:auto}h1{font-size:20px;font-weight:normal;margin:0 0 16px}
select,input{font:inherit;padding:6px;margin:0 12px 8px 0;max-width:100%}
img{display:block;max-width:100%;height:auto}img[hidden]{display:none}
</style><main><h1>RNA indels by target</h1>
<label>Dataset <select id="cohort"><option value="company">Company merged</option><option value="chip01">Chip 1</option></select></label>
<label>Find <input id="search" type="search" placeholder="Target name or ID"></label>
<label>Target <select id="target"></select></label>
<p id="count"></p><a id="full" title="Open full-resolution PNG"><img id="figure" alt="PCA of cells targeting the selected ORF"></a>
</main><script>
const entries = __ENTRIES__;
const cohort = document.querySelector('#cohort');
const search = document.querySelector('#search');
const target = document.querySelector('#target');
const figure = document.querySelector('#figure');
const full = document.querySelector('#full');
const count = document.querySelector('#count');
let shown = [];
function display() {
    const item = shown.find(row => row.id === target.value);
    figure.hidden = !item;
    if (!item) { figure.removeAttribute('src'); full.removeAttribute('href'); return; }
    figure.src = item.png;
    full.href = item.png;
}
function filter(prefer = target.value) {
    const query = search.value.toLowerCase().trim();
    shown = entries.filter(row => row.dataset === cohort.value && (row.name + ' ' + row.id).toLowerCase().includes(query));
    target.replaceChildren(...shown.map(row => new Option(row.name + ' (' + row.indels + ' indels)', row.id)));
    if (shown.some(row => row.id === prefer)) target.value = prefer;
    count.textContent = shown.length + ' targets';
    display();
}
cohort.onchange = () => filter();
search.oninput = () => filter();
target.onchange = display;
filter('PB.26158.49_chr2_1_2485_32_236_uORF_ATG');
</script></html>'''.replace('__ENTRIES__', data), encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('coordinates', 'cases', 'labels', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--commit', required=True)
    parser.add_argument('--preview-target', action='append', default=[])
    args = parser.parse_args()
    with args.cases.open(encoding='utf-8', newline='') as handle:
        cases = {r['dataset'] + '/' + r['target_gene_id']: r for r in csv.DictReader(handle)}
    labels = json.loads(args.labels.read_text(encoding='utf-8'))
    variants = defaultdict(list)
    for dataset in COHORT:
        for event in labels[dataset]['variants']:
            if event['kind'] in ('INS', 'DEL') and event['guide_alt'] > 0:
                key = dataset + '/' + event['target_gene_id']
                assert key in cases
                variants[key].append(event)

    font = Path('C:/Windows/Fonts/arial.ttf')
    configure(font if font.is_file() else None)
    entries = []
    counts = Counter()
    with zipfile.ZipFile(args.coordinates) as archive:
        sources = json.loads(archive.read('manifest.json'))['files']
        assert set(cases) == set(sources)
        args.output.mkdir(parents=True, exist_ok=False)
        for key, row in sorted(cases.items()):
            if args.preview_target and row['target_gene_id'] not in args.preview_target:
                continue
            npz = archive.read(key + '.npz')
            assert hashlib.sha256(npz).hexdigest() == sources[key]['subset_sha256']
            events = sorted(variants[key], key=lambda v: (v['guide'], v['pos0'], v['ref'], v['alt']))
            filename = key + '.png'
            destination = args.output / filename
            destination.parent.mkdir(exist_ok=True)
            with np.load(io.BytesIO(npz), allow_pickle=False) as arrays:
                annotations = draw_case(row, arrays, events, destination)
            assert len(annotations) == len(events)
            entries.append({'dataset': row['dataset'], 'id': row['target_gene_id'], 'name': target_display(row),
                            'records': int(row['RNA_records']), 'indels': len(events), 'png': filename,
                            'sha256': hashlib.sha256(destination.read_bytes()).hexdigest(),
                            'source_PCA_sha256': sources[key]['source_sha256'], 'changes': annotations})
            counts[row['dataset']] += 1
            if counts[row['dataset']] % 25 == 0:
                print(json.dumps({'figures_written': dict(counts)}), flush=True)
    assert entries
    if not args.preview_target:
        assert len(entries) == len(cases)
        assert sum(r['indels'] for r in entries) == sum(map(len, variants.values()))
    gallery(args.output / 'index.html', entries)
    receipt = {'code_commit': args.commit, 'complete': not args.preview_target,
               'inputs_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in (args.coordinates, args.cases, args.labels)},
               'PCA_scope': 'all RNA records for the same target', 'background_NTC_records': 0,
               'figures': dict(counts), 'indels': sum(r['indels'] for r in entries), 'cases': entries}
    (args.output / 'render_receipt.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({key: receipt[key] for key in ('complete', 'figures', 'indels')}), flush=True)


if __name__ == '__main__':
    main()
