import sys
import os
import subprocess


controllerNode = sys.argv[1]
controllerPort = sys.argv[2]
BASEPATH  = sys.argv[3]
INDEXDIRPATH = sys.argv[4]

# rsync code
args = ["rsync", '-a', INDEXDIRPATH, '/var/tmp/']
_ = subprocess.call(args)

returnCode = 1
executable = BASEPATH + '/mrsfast2/mrsfast'

while (returnCode == 1):
	args = [executable, '-N', controllerNode, '-P', controllerPort]
	returnCode = subprocess.call(args)

