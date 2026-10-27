"""Apply a hash-gated log guard adapter to an unused Martian source copy."""
from hashlib import sha256
from pathlib import Path
import argparse


def adapt(raw):
    expected = Path(__file__).with_name('log_source_sha256.txt').read_text().strip()
    if sha256(raw).hexdigest() != expected:
        raise ValueError('Pinned upstream log.go differs')
    text = raw.decode()
    text = text.replace('\t"os"\n', '\t"os"\n\t"sync"\n', 1)
    text = text.replace('type Logger struct {\n', 'type Logger struct {\n\tfileMutex sync.Mutex\n', 1)
    for signature in ['func (logger *Logger) Write(msg []byte) (int, error) {',
                      'func (logger *Logger) WriteString(msg string) (int, error) {']:
        assert text.count(signature) == 1
        text = text.replace(signature, signature+'\n\tlogger.fileMutex.Lock()\n\tdefer logger.fileMutex.Unlock()', 1)
    anchor = '\tw := LOGGER.fileWriter\n'
    assert text.count(anchor) == 1
    text = text.replace(anchor, '\tLOGGER.fileMutex.Lock()\n\tdefer LOGGER.fileMutex.Unlock()\n'+anchor, 1)
    anchor = '\t} else if !os.SameFile(info, currentInfo) {\n'
    assert text.count(anchor) == 1
    text = text.replace(anchor, anchor+'\t\tif proofErr := verifyCifsLogIdentity(f, info, currentInfo); proofErr == nil {\n\t\t\treturn nil\n\t\t}\n', 1)
    return text.encode()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('unused_source', type=Path)
    args = parser.parse_args()
    source = args.unused_source.resolve(strict=True)
    if (source/'bin/mrp').exists():
        raise ValueError('Use a fresh source copy without a runtime; never modify an active view')
    target = source/'martian/util/log.go'
    target.write_bytes(adapt(target.read_bytes()))
    for name in ['cifs_log_identity_linux.go', 'cifs_log_identity_other.go', 'cifs_log_identity_linux_test.go']:
        out = target.parent/name
        if out.exists():
            raise ValueError('Adapter target already exists')
        out.write_bytes(Path(__file__).with_name(name).read_bytes())


if __name__ == '__main__':
    main()
