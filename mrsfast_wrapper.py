import sys
import os
import signal
import subprocess
import argparse

global pid

def SIGINT_handler(s, frame):
    os.kill(pid, signal.SIGINT)
    print 'Mapper EXITING with return code 0!'
    sys.exit(1)

if __name__ == '__main__':
    signal.signal(signal.SIGINT, SIGINT_handler) # handle Ctrl-C gracefully
    parser = argparse.ArgumentParser()
    parser.add_argument("controllerNode")
    parser.add_argument("controllerPort")
    parser.add_argument("BASEPATH")
    parser.add_argument("--index-dir", required=False, default=None)
    parser.add_argument("--temp", required=False, default='/var/tmp/')
    args = parser.parse_args()
    # rsync code
    print args
    if args.index_dir:
        _ = subprocess.call(["rsync", '-a', args.index_dir, args.temp])

    returnCode = 1
    executable = args.BASEPATH + '/mrsfast2/mrsfast'

    exec_args = [executable, '-N', args.controllerNode, '-P', args.controllerPort]

    while (returnCode == 1):
        #returnCode = subprocess.call(exec_args)
        p = subprocess.Popen(exec_args)
        pid = p.pid
        returnCode = p.wait()

