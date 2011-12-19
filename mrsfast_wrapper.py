import sys
import os
import subprocess


controllerNode = sys.argv[1]
controllerPort = sys.argv[2]
BASEPATH  = sys.argv[3]

returnCode = 1
executable = BASEPATH + '/mrsfast/mrsfast'

while (returnCode == 1):
	args = [executable, '-N', controllerNode, '-P', controllerPort]
	returnCode = subprocess.call(args)

