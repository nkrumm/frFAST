import sys
import os
import subprocess
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("controllerNode")
    parser.add_argument("controllerPort")
    parser.add_argument("BASEPATH")
    parser.add_argument("INDEXDIRPATH", required=False, default=None)
    args = parser.parse_args()
    # rsync code

    if args.INDEXDIRPATH:
        args = ["rsync", '-a', args.INDEXDIRPATH, '/var/tmp/']
        _ = subprocess.call(args)

    returnCode = 1
    executable = args.BASEPATH + '/mrsfast2/mrsfast'

    args = [executable, '-N', args.controllerNode, '-P', args.controllerPort]

    while (returnCode == 1):    
        returnCode = subprocess.call(args)
