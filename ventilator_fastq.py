import sys
import zmq
import numpy as np
import os
#import pysam
from subprocess import Popen, PIPE

import time

# def gziplines(fname): 
# 	f = Popen(['zcat', fname], stdout=PIPE)
# 	for line in f.stdout:
# 		yield line

context = zmq.Context()

# connect to controller

controlleraddress = sys.argv[1]
print "Connecting to controller @ ", controlleraddress

requester = context.socket(zmq.REQ);
requester.connect("tcp://" + controlleraddress + ":5555")

nodeaddress = os.uname()[1]

requester.send("VENT REGISTER " + nodeaddress)
_ = requester.recv()



# get bam file location

requester.send("VENT GETSOURCE")
source_filename = requester.recv()

print "got source: " , source_filename


# get bam file location

requester.send("VENT GETNUMREADS")
maximum_reads_per_destination = int(requester.recv())


# WAIT for mappers to start and "ventilate" signal from controller
requester.send("VENT READY")
_ = requester.recv()



f_src=open(source_filename)
total_reads = 0
read_counter = 0
endOfFile = False 

while not endOfFile:
	requester.send("VENT GETDEST")
	node = requester.recv()
	while node == "0":
		requester.send("VENT GETDEST")
		node = requester.recv()
		print "waiting for mapper to become available ..."
		time.sleep(2)
	
	print "received as mapper node: ", node
	
	readstreamer = context.socket(zmq.PUSH)
	readstreamer.bind("tcp://*:8000")
	readstreamer.send(str(maximum_reads_per_destination))
	read_counter = 0
	total_reads = 0
	
	while read_counter < maximum_reads_per_destination:
		readname = f_src.readline().strip("\n")
		
		if not readname:
			endOfFile = True
			break
		
		seq = f_src.readline().strip("\n")
		msg = readname + " " + seq		
		readstreamer.send(msg)
		read_counter += 1
		total_reads += 1	
		#print total_reads, read_counter, msg
	
	readstreamer.send("DONE 0")
	print "Done ventilating to " + node
	requester.send("VENT DONE " + str(read_counter) + " " + str(total_reads) + " " + node)
	_ = requester.recv()
	readstreamer.close()
	time.sleep(3);


print "Finished all ventilating all reads!"
requester.send("VENT FINISHED " + str(read_counter) + " " + str(total_reads))
_ = requester.recv()
f_src.close()