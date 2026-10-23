import time

def main(args, outs):
    with open(args.started, 'w') as f:
        f.write('waiting\n')
    time.sleep(20)
    with open(outs.result, 'w') as f:
        f.write(str(args.value + 1) + '\n')
