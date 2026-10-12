"""Execute the versioned QC notebook into a user-selected result directory."""
from __future__ import annotations

import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import subprocess
import sys

import nbformat
from nbclient import NotebookClient


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    config_path = args.config.resolve()
    config = json.loads(config_path.read_text(encoding='utf-8'))
    output_root = Path(config['output_root']).resolve()
    if output_root.is_relative_to(repo):
        raise ValueError('Analysis results must be outside the code repository')
    output = output_root / 'provenance/01_lab_chip01_analysis.executed.ipynb'
    output.parent.mkdir(parents=True, exist_ok=True)
    lock = output.parent / 'notebook_execution.lock'
    with lock.open('x', encoding='utf-8') as handle:
        handle.write(json.dumps({'pid': os.getpid(), 'started': datetime.now().isoformat()}))
    notebook = nbformat.read(repo / 'notebooks/01_lab_chip01_analysis.ipynb', as_version=4)
    try:
        subprocess.run([sys.executable, '-m', 'ipykernel', 'install', '--prefix', sys.prefix,
                        '--name', 'perturbseq-analysis', '--display-name', 'Perturb-seq analysis'], check=True)
        os.environ['JUPYTER_PATH'] = str(Path(sys.prefix) / 'share/jupyter')
        client = NotebookClient(notebook, timeout=7200, kernel_name='perturbseq-analysis',
                                resources={'metadata': {'path': str(repo)}})
        environment = dict(os.environ, PERTURBSEQ_ANALYSIS_CONFIG=str(config_path),
                           OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2')
        client.execute(env=environment)
        temporary = output.with_suffix('.tmp.ipynb')
        nbformat.write(notebook, temporary)
        temporary.replace(output)
        print(f'NOTEBOOK_COMPLETE {output}', flush=True)
    except BaseException:
        nbformat.write(notebook, output.with_name('01_lab_chip01_analysis.failed.ipynb'))
        raise
    finally:
        lock.unlink(missing_ok=True)


if __name__ == '__main__':
    main()
