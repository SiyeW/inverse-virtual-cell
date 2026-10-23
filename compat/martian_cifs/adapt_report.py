"""Build a local-only reporting overlay without redistributing vendor source."""
import ast,hashlib

REPORT_PATH='mro/rna/stages/counter/summarize_reports/__init__.py'

def adapt(source,expected_sha256):
    if hashlib.sha256(source).hexdigest()!=expected_sha256:raise ValueError('Reporting source version differs')
    text=source.decode('utf8');tree=ast.parse(text)
    def field(node,owner):
        return isinstance(node,ast.Attribute) and isinstance(node.value,ast.Name) and node.value.id==owner and node.attr=='feature_reference'
    matches=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)
             and isinstance(n.func.value,ast.Name) and n.func.value.id=='shutil' and n.func.attr=='copy'
             and len(n.args)==2 and not n.keywords and field(n.args[0],'args') and field(n.args[1],'outs')]
    if len(matches)!=1:raise ValueError('Expected exactly one report attachment copy')
    old=ast.get_source_segment(text,matches[0])
    if text.count(old)!=1 or text.count('import shutil\n')!=1:raise ValueError('Unexpected report import/call layout')
    new='cifs_content_copy.checked_copy('+ast.get_source_segment(text,matches[0].args[0])+', '+ast.get_source_segment(text,matches[0].args[1])+')'
    result=text.replace(old,new,1).replace('import shutil\n','import shutil\nimport cifs_content_copy\n',1)
    compile(result,'internal_reporting_overlay','exec')
    return result.encode('utf8')
