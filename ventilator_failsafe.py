import sys
import zmq
import numpy as np
import os
import pysam
from subprocess import Popen, PIPE
import string
import glob
import time

def errorMsg(msg):
	requester.send(msg)
	_ = requester.recv()
	exit(0)

context = zmq.Context(8)

# connect to controller

controlleraddress = sys.argv[1]
controllerport = int(sys.argv[2])

print "Connecting to controller @ " +str(controlleraddress) + ":" + str(controllerport)

VENT2MAPPER_PORT = controllerport +2 #formerly  :8000

print "Will be using port %d as vent-->mapper port (formerly :8000)" % VENT2MAPPER_PORT

try:
	requester = context.socket(zmq.REQ);
	requester.connect("tcp://" + controlleraddress + ":" +str(controllerport))
except zmq.core.error.ZMQError:
	# port already in use!
	# obviouslly cannot communicate via this port, must just exit.
	print "ERROR VENT controller port not available!"
	sys.exit(0)

nodeaddress = os.uname()[1]

requester.send("VENT REGISTER " + nodeaddress)
_ = requester.recv()



def gzipfolder_line(folder_path):
	f = Popen('zcat ' + folder_path, shell=True, stdout=PIPE)
	read = pysam.AlignedRead()
	for line in f.stdout:
		read.qname = line.strip("\n")
		read.seq = f.stdout.next().strip("\n")
		_ = f.stdout.next()
		_ = f.stdout.next()
		yield read

# get bam file location

requester.send("VENT GETSOURCE")
source_filename = requester.recv()

print "got source: " , source_filename

## attempt to figure out file type:
if source_filename[-3:].lower() == 'bam':
	source_type = 'bam'
	try:
		s = pysam.Samfile(source_filename,'rb')
	except IOError:
		errorMsg("ERROR VENT cannot open BAM file: " + source_filename)
	
	source = s.fetch(until_eof=True)
	
elif os.path.isdir(source_filename):
	source_type = 'fastq_folder'
	source_filename = source_filename + "/*.gz"
	num_files = len(glob.glob(source_filename))
	if num_files == 0:
		errorMsg("ERROR VENT no FASTQ files found in folder (or folder does not exist): " + source_filename)
	try:
		source = gzipfolder_line(source_filename)
	except:
		errorMsg("ERROR VENT cannot open fastq.gz folder: " + source_filename)
	
else:
	errorMsg("ERROR VENT cannot determine file type: " + source_filename)


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
	rc_offset_start_36 = 0
	end_split = 72
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
start_positions  = {}

current_chunkID = 1

while not endOfFile:
	requester.send("VENT GETDEST")
	msg = requester.recv()
	while msg == "0":
		requester.send("VENT GETDEST")
		msg = requester.recv()
		print "waiting for mapper to become available ..."
		time.sleep(2)
	
	node,mapperID,chunkID = msg.split(" ") # TODO multithread need to receive port, as well as BAM file offset to use.
	
	mapperID = int(mapperID)
	chunkID = int(chunkID)
	
	print "Received node: %s, mapperID: %d, chunkID %d" % (node, mapperID, chunkID)
	
	try:
		readstreamer = context.socket(zmq.PUSH)
		readstreamer.bind("tcp://*:" + str(VENT2MAPPER_PORT))
	except zmq.core.error.ZMQError:
		# port already in use!
		errorMsg("ERROR VENT mapper readstream port not available!")
	
	
	readstreamer.send(str(maximum_reads_per_destination))
	read_counter = 0
	#total_reads = 0
	t1 = time.time()
		
	if current_chunkID == chunkID: # we expect an incremented count for chunkID
		if source_type == 'bam':
			start_positions[chunkID] = s.tell()
		elif source_type == 'fastq_folder':
			pass # todo
	else: # however, if the chunkID does not increment as expected, we seek back to the requested chunkID
		s.seek(start_positions[chunkID])
		print "Retrying... chunkID = %d" % chunkID
	
	while read_counter < maximum_reads_per_destination:
		try:
			read = source.next()
		except StopIteration:
			endOfFile = True
			break
		
		if read.is_reverse:
			if read.rlen == 36:
				readstreamer.send("%s %s" % (read.qname, reverseComplement(read.seq[rc_offset_start_36:])))
				read_counter += 1
			elif read.rlen == 50:
				readstreamer.send("%s %s" % (read.qname, reverseComplement(read.seq[rc_offset_start_50:])))
				read_counter += 1				
			elif read.rlen == 76:
				seq = reverseComplement(read.seq[rc_offset_start_76:])
				readstreamer.send("%s.%s %s" % (read.qname, "01", seq[0:split_length]))
				readstreamer.send("%s.%s %s" % (read.qname, "02", seq[split_length:]))	
				read_counter += 2
			else:
				# send error message
				requester.send("WARNING VENT found non-standard read length: "  + str(read.rlen) + " -- " + str(read))
				_ = requester.recv()
				#exit(0)
		else:
			if read.rlen == 36:
				readstreamer.send("%s %s" % (read.qname,read.seq[0:split_length]))
				read_counter += 1
			elif read.rlen == 50:
				readstreamer.send("%s %s" % (read.qname,read.seq[0:split_length]))
				read_counter += 1
			elif read.rlen == 76:
				readstreamer.send("%s.%s %s" % (read.qname,"01",read.seq[0:split_length]))
				readstreamer.send("%s.%s %s" % (read.qname,"02",read.seq[split_length:end_split]))
				read_counter += 2
			else:
				# send error message
				requester.send("WARNING VENT found non-standard read length: " + str(read.rlen) + " -- " + str(read))
				_ = requester.recv()
				#exit(0)
	
	
	current_chunkID +=1	
		
	print "Done ventilating to " + node + " in " + str(time.time() -t1) + " s"
	#sys.stdout.flush()
	
	readstreamer.send("DONE 0")
	readstreamer.close()
	
	requester.send("VENT DONE " + str(read_counter) + " " + str(read_counter) + " " + node + " " + str(mapperID))
	_ = requester.recv()
	
	time.sleep(1);


print "Finished all ventilating all reads!"
requester.send("VENT FINISHED " + str(read_counter) + " " + str(read_counter))
_ = requester.recv()