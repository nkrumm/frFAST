import sys
import zmq
import numpy as np
import os
#import pysam
from subprocess import Popen, PIPE
import string

import time

context = zmq.Context(8)

# connect to controller

controlleraddress = sys.argv[1]
controllerport = int(sys.argv[2])

print "Connecting to controller @ " +str(controlleraddress) + ":" + str(controllerport)

VENT2MAPPER_PORT = controllerport +2 #formerly  :8000

print "Will be using port %d as vent-->mapper port (formerly :8000)" % VENT2MAPPER_PORT

requester = context.socket(zmq.REQ);
requester.connect("tcp://" + controlleraddress + ":" +str(controllerport))

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

def samlines(fname): 
	f = Popen(['samtools', 'view', fname], stdout=PIPE)
	while 1:
		a = f.stdout.readlines(10000000) # 1000000 = 2400 lines, 10000000 = 24000
		for i in range(len(a)):
			yield a[i]

samfile = samline(source_filename)


# WAIT for mappers to start and "ventilate" signal from controller
requester.send("VENT READY")
_ = requester.recv()

complement = string.maketrans('atcgnATCGN', 'tagcnTAGCN')

def reverseComplement(sequence):
	return sequence.translate(complement)[::-1]


total_reads = 0
read_counter = 0
endOfFile = False 

rc_flag = 16 #check for reverse complement in BAM file
unmapped_flag = 141 # this is the second read of the unmapped reads!

while not endOfFile:
	requester.send("VENT GETDEST")
	msg = requester.recv()
	while msg == "0":
		requester.send("VENT GETDEST")
		msg = requester.recv()
		print "waiting for mapper to become available ..."
		time.sleep(2)
	
	node,mapperID = msg.split(" ") # TODO multithread need to receive port, as well as BAM file offset to use.
	
	# TODO, multithread here. will need to pull in additional file pipes to read BAM file.
	
	print "received as mapper node: ", node
	
	readstreamer = context.socket(zmq.PUSH)
	readstreamer.bind("tcp://*:" + str(VENT2MAPPER_PORT)) # TODO, set up the port based on received information
	readstreamer.send(str(maximum_reads_per_destination))
	read_counter = 0
	total_reads = 0
	t1 = time.time()
	
	while read_counter < maximum_reads_per_destination:
		try:
			read = samfile.next()
		except StopIteration:
			endOfFile = True
			break
		
		read = read.split("\t",10)
		
		if (int(read[1]) & rc_flag):# needs to be reverse complemented.
			seq = reverseComplement(read[9])
		else:
			seq = read[9]
		
		if number_splits == 1:
			msg = "%s %s" % (read[0],seq[0:split_length])
			readstreamer.send(msg)
			read_counter += 1
			total_reads += 1	
			#print total_reads, read_counter, msg
		elif number_splits == 2:
			#SPLIT 1:
			msg = "%s.1 %s" % (read[0],seq[0:split_length])
			readstreamer.send(msg)
			#SPLIT 2:
			msg = "%s.2 %s" % (read[0],seq[split_length:(split_length*2)])
			readstreamer.send(msg)
			
			read_counter += 2
			total_reads += 2	
			
			
	readstreamer.send("DONE 0")
	print "Done ventilating to " + node + " in " + str(time.time() -t1) + " s"
	requester.send("VENT DONE " + str(read_counter) + " " + str(total_reads) + " " + node + " " + str(mapperID))
	_ = requester.recv()
	readstreamer.close()
	time.sleep(0.1);


print "Finished all ventilating all reads!"
requester.send("VENT FINISHED " + str(read_counter) + " " + str(total_reads))
_ = requester.recv()