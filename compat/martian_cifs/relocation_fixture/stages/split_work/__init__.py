import time
from pathlib import Path

def split(args):
    return {'chunks':[{'index':0},{'index':1}]}

def main(args,outs):
    if args.index==1:
        with Path(args.marker).open('x') as f:f.write('completed chunk ran once\n')
    else:
        while not Path(args.release).exists():time.sleep(.2)
    Path(outs.piece).write_text(str(args.index+20)+'\n')

def join(args,outs,chunk_defs,chunk_outs):
    value=sum(int(Path(x.piece).read_text()) for x in chunk_outs)
    Path(outs.result).write_text(str(value)+'\n')
