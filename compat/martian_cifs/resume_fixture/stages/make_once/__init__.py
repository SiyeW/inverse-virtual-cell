def main(args, outs):
    # Re-execution is an error: successful resume must reuse this completed step.
    with open(args.marker, 'x') as f:
        f.write('executed once\n')
    outs.value = 41
