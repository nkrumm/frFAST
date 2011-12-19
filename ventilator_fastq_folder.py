
import sys
import zmq
import numpy as np
import os
#import pysam
from subprocess import Popen, PIPE

import time

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


# get number reads per mapper

requester.send("VENT GETNUMREADS")
maximum_reads_per_destination = int(requester.recv())


requester.send("VENT GETSPLITPARAMS")
data = requester.recv()
split_length, read_length = data.split(" ")
split_length = int(split_length)
read_length = int(read_length)

number_splits = read_length / split_length

def samline(fname): 
	f = Popen(['samtools', 'view', fname], stdout=PIPE)
	for line in f.stdout:
		yield line

def gzipfolder_line(folder_path):	
	folder_path = folder_path + "*.gz"
	f = Popen('zcat ' + folder_path, shell=True, stdout=PIPE)
	for line in f.stdout:
		yield line



samfile = gzipfolder_line(source_filename)


# WAIT for mappers to start and "ventilate" signal from controller
requester.send("VENT READY")
_ = requester.recv()


total_reads = 0
read_counter = 0
endOfFile = False 

while not endOfFile:
	requester.send("VENT GETDEST")
	msg = requester.recv()
	while msg == "0":
		requester.send("VENT GETDEST")
		msg = requester.recv()
		print "waiting for mapper to become available ..."
		time.sleep(2)
	
	node,mapperID = msg.split(" ")
	
	print "received as mapper node: ", node
	
	readstreamer = context.socket(zmq.PUSH)
	readstreamer.bind("tcp://*:8000")
	readstreamer.send(str(maximum_reads_per_destination))
	read_counter = 0
	total_reads = 0
	
	while read_counter < maximum_reads_per_destination:
		try:
			readname = samfile.next().strip("\n")
			seq = samfile.next().strip("\n")
			_ = samfile.next() # ignore next two lines
			_ = samfile.next()
		except StopIteration:
			endOfFile = True
			break
		
# 		read = read.split("\t",10)
# 		readname = read[0]
# 		seq = read[9]
		if number_splits == 1:
			msg = readname + " " + seq		
			readstreamer.send(msg)
			read_counter += 1
			total_reads += 1	
			#print total_reads, read_counter, msg
		elif number_splits == 2:
			#SPLIT 1:
			msg = readname + ".1 " + seq[0:split_length]
			readstreamer.send(msg)
			#SPLIT 1:
			msg = readname + ".2 " + seq[split_length:(split_length*2)]
			readstreamer.send(msg)
			
			read_counter += 2
			total_reads += 2	
			
			
	readstreamer.send("DONE 0")
	print "Done ventilating to " + node
	requester.send("VENT DONE " + str(read_counter) + " " + str(total_reads) + " " + node + " " + str(mapperID))
	_ = requester.recv()
	readstreamer.close()
	time.sleep(1);


print "Finished all ventilating all reads!"
requester.send("VENT FINISHED " + str(read_counter) + " " + str(total_reads))
_ = requester.recv()