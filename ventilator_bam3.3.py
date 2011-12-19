import sys
import zmq
import numpy as np
import os
import pysam
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

try:
	s = pysam.Samfile(source_filename,'rb')
except IOError:
	requester.send("ERROR VENT cannot open BAM file: " + source_filename)
	_ = requester.recv()
	exit(0)

samfile = s.fetch(until_eof=True)

# get number reads per mapper
requester.send("VENT GETNUMREADS")
maximum_reads_per_destination = int(requester.recv())


requester.send("VENT GETSPLITPARAMS")
data = requester.recv()
split_length, read_length = data.split(" ")
split_length = int(split_length)
if read_length == 'AUTO':
	rc_offset_start_50 = 14
	rc_offset_start_76 = 4
	end_split = 72
	#read = samfile.next()
	#read_length = len(read.seq)
	#samfile = s.fetch(until_eof=True) # rewind iterator
else:
	read_length = int(read_length)


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
	#total_reads = 0
	t1 = time.time()
	
	while read_counter < maximum_reads_per_destination:
		try:
			read = samfile.next()
		except StopIteration:
			endOfFile = True
			break
		
		if read.is_reverse:
			if read.rlen == 50:
				readstreamer.send("%s %s" % (read.qname, reverseComplement(read.seq[rc_offset_start_50:])))
				read_counter += 1				
			elif read.rlen == 76:
				seq = reverseComplement(read.seq[rc_offset_start_76:])
				readstreamer.send("%s.%s %s" % (read.qname, "01", seq[0:split_length]))
				readstreamer.send("%s.%s %s" % (read.qname, "02", seq[split_length:]))	
				read_counter += 2
			else:
				# send error message
				requester.send("ERROR VENT found non-standard read length: "  + str(read.rlen) + " -- " + str(read))
				_ = requester.recv()
				exit(0)
		else:
			if read.rlen == 50:
				readstreamer.send("%s %s" % (read.qname,read.seq[0:split_length]))
				read_counter += 1
			elif read.rlen == 76:
				readstreamer.send("%s.%s %s" % (read.qname,"01",read.seq[0:split_length]))
				readstreamer.send("%s.%s %s" % (read.qname,"02",read.seq[split_length:end_split]))
				read_counter += 2
			else:
				# send error message
				requester.send("ERROR VENT found non-standard read length: " + str(read.rlen) + " -- " + str(read))
				_ = requester.recv()
				exit(0)
	
	readstreamer.send("DONE 0")
	print "Done ventilating to " + node + " in " + str(time.time() -t1) + " s"
	requester.send("VENT DONE " + str(read_counter) + " " + str(read_counter) + " " + node + " " + str(mapperID))
	_ = requester.recv()
	readstreamer.close()
	time.sleep(1);


print "Finished all ventilating all reads!"
requester.send("VENT FINISHED " + str(read_counter) + " " + str(read_counter))
_ = requester.recv()