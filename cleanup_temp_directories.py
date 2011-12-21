# execute_on_all_nodes.py
# MUST RUN ON EEEK
# author: Nik Krumm
# 
import time
import sys
import os
import subprocess
from xml.dom.minidom import parse, parseString

# Queue on which to execute the commands. Commands in commandList will be executed on each of the queue's nodes!
target_queues = ["all.q","prod.q"]

# these are the commands, executed in order!
commandList = ['for f in `ls /var/tmp/*.h5`; do if test `find "$f" -mmin +120`; then rm -f $f; fi; done']


cmd = "qhost -q -xml"
out,e = subprocess.Popen(cmd.split(" "),stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()

dom = parseString( out )

target_queue_hosts = []

for host in dom.getElementsByTagName("host"):
	hostname = host.getAttribute("name")
	for q in host.getElementsByTagName("queue"):
		if str(q.getAttribute("name")) in target_queues:
			target_queue_hosts.append(hostname)

for hostname in target_queue_hosts:
	t1 = time.time()
	for cmd in commandList:
		rshcmd = "rsh " + hostname + " " + cmd
		out,e = subprocess.Popen(rshcmd.split(" "),stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
		print hostname + "> " + cmd
		print out, e
	print time.time()-t1