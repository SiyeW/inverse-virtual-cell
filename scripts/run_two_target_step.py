"""Run one pinned, resource-bounded stage of the two-target RNA pilot."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import socket
import subprocess
import time


def startticks(pid):
    value = Path(f'/proc/{pid}/stat').read_text()
    return int(value[value.rfind(')') + 2:].split()[19])


def atomic(path, value):
    partial = path.with_suffix('.receiving')
    partial.write_text(json.dumps(value, indent=2) + '\n')
    partial.replace(path)


def run(config):
    if os.getuid() != 2601 or socket.gethostname() != 'rhesusbase':
        raise ValueError('This stage belongs to the adopted laboratory host')
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[:2])
    native = Path(config['native'])
    if not native.is_dir() or (native / 'run.exit').exists():
        raise ValueError('Expected a new, prepared stage directory')
    lock = native / 'active.lock'
    descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    owner = dict(pid=os.getpid(), uid=os.getuid(), startticks=startticks(os.getpid()),
                 stage=config['stage'], commit=config['commit'])
    os.write(descriptor, json.dumps(owner).encode())
    os.close(descriptor)
    child = None
    begin = time.monotonic()

    def state(stage, **details):
        atomic(native / 'status.json', dict(stage=stage, at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),
               elapsed_seconds=round(time.monotonic() - begin, 1), **details))

    try:
        state('verify_pinned_stage')
        for name, expected in config['sha256'].items():
            with Path(name).open('rb') as stream:
                actual = hashlib.file_digest(stream, 'sha256').hexdigest()
            if actual != expected:
                raise ValueError('Pinned source differs: ' + name)
        environment = os.environ | config.get('environment', {})
        environment.update(OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2',
                           NUMEXPR_NUM_THREADS='2', NUMBA_NUM_THREADS='2')

        def memory_limit():
            maximum = 32 * 1024**3
            resource.setrlimit(resource.RLIMIT_AS, (maximum, maximum))

        with (native / 'stdout.log').open('w') as stdout, (native / 'stderr.log').open('w') as stderr:
            child = subprocess.Popen(config['command'], env=environment, stdout=stdout, stderr=stderr,
                                     preexec_fn=memory_limit)
            identity = dict(pid=child.pid, uid=os.getuid(), startticks=startticks(child.pid))
            while child.poll() is None:
                phase_file = Path(config['phase_file'])
                phase = phase_file.read_text().strip() if phase_file.is_file() else 'starting'
                state('running', phase=phase, owner=identity, command=config['command'])
                if time.monotonic() - begin > config.get('protection_seconds', 43200):
                    raise TimeoutError('Stage protection reached; completed files are retained')
                time.sleep(10)
        if child.returncode:
            state('failed_pending_recovery', exit=child.returncode)
            return child.returncode
        for name in config['required_results']:
            if not Path(name).is_file():
                raise ValueError('Stage did not publish required result: ' + name)
        state('complete_pending_acceptance', exit=0)
        return 0
    except BaseException as error:
        if child is not None and child.poll() is None:
            child.terminate()
            try:
                child.wait(timeout=20)
            except subprocess.TimeoutExpired:
                child.kill()
                child.wait(timeout=20)
        state('failed_pending_recovery', error=f'{type(error).__name__}: {error}')
        raise
    finally:
        lock.unlink(missing_ok=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, required=True)
    arguments = parser.parse_args()
    raise SystemExit(run(json.loads(arguments.config.read_text())))
