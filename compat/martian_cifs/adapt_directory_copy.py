"""Hash-gated, one-call adaptation of a local vendor cr_io.py copy."""
import ast
import hashlib

CR_IO_PATH = 'lib/python/cellranger/cr_io.py'


def adapt(source, expected_sha256):
    if hashlib.sha256(source).hexdigest() != expected_sha256:
        raise ValueError('Unexpected cr_io.py source identity')
    tree = ast.parse(source)
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef)
             and n.name == 'hardlink_with_fallback']
    calls = [n for n in ast.walk(nodes[0]) if isinstance(n, ast.Call)
             and isinstance(n.func, ast.Attribute) and isinstance(n.func.value, ast.Name)
             and n.func.value.id == 'shutil' and n.func.attr == 'copytree'] if len(nodes) == 1 else []
    text = source.decode()
    old = 'shutil.copytree(src, dst, copy_function=_hardlink_file_with_fallback, dirs_exist_ok=True)'
    if len(calls) != 1 or text.count(old) != 1 or text.count('import shutil\n') != 1:
        raise ValueError('Unexpected directory publication call layout')
    result = text.replace('import shutil\n', 'import shutil\nimport cifs_directory_copy\n', 1)
    result = result.replace(old, old.replace('shutil.copytree', 'cifs_directory_copy.copytree'), 1)
    ast.parse(result)
    return result.encode()
