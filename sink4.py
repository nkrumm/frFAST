# required modules: 
# module import zlib hdf5/1.8.4 
# module import hdf5/1.8.4 (or 1.8.2 on some nodes)

import numpy as np
import math
import operator
import random
import time
import sys
import itertools
import re
import os
from tables import *
from collections import defaultdict
from subprocess import Popen, PIPE
import zmq

context = zmq.Context(8)


controlleraddress = sys.argv[1]
controllerport = int(sys.argv[2])

print "Connecting to controller @ " +str(controlleraddress) + ":" + str(controllerport)

MAPPER2SINK_PORT = controllerport +3 #formerly  :6000

print "Will be using port %d as mapper-->sink port (formerly :6000)" % MAPPER2SINK_PORT


requester = context.socket(zmq.REQ);
requester.connect("tcp://" + controlleraddress + ":" + str(controllerport))

nodeaddress = os.uname()[1]

requester.send("SINK REGISTER " + nodeaddress)
_ = requester.recv()

insock = context.socket(zmq.PULL)
insock.bind("tcp://*:"+ str(MAPPER2SINK_PORT))

quit_signal = False

total_mappings = 0

requester.send("SINK GETDEST " + nodeaddress)
outfilename = requester.recv()
print "outfilename: " , outfilename

requester.send("SINK READY")
_ = requester.recv()


##### SET UP OUTARRAYS
MAX_READS = 2000000

mapped_chrs = np.zeros([MAX_READS],dtype=np.uint8)
mapped_pos = np.zeros([MAX_READS],dtype=np.uint32)
mapped_perfect = np.zeros([MAX_READS],dtype=np.bool)
mapped_mapperID = np.zeros([MAX_READS],dtype=np.uint8)
counter = 0

#### SET UP TEMP HDF5 
temp_hdf5_filename = "/var/tmp/frFAST.tmp.h5"
h5file = openFile(temp_hdf5_filename, mode = "w", title = "Temp")
group = h5file.createGroup("/",'mappings','Mappings')

class mapped_read(IsDescription):
	pos = UInt32Col(pos=0)
	perfect_match = BoolCol(pos=1)

tables = []
for chr in range(1,25):
	tables.append(h5file.createTable(group,"chr"+str(chr),mapped_read,"chr"+str(chr)))#,expectedrows=mapping_exprows[chr]))

for chr in range(31,55):
	tables.append(h5file.createTable(group,"chr"+str(chr),mapped_read,"chr"+str(chr)))#,expectedrows=mapping_exprows[chr]))

h5file.flush()

while 1: 
	data_msg = insock.recv()
	try:
		mapperID, chr, position, edit = data_msg.split(" ")
	except ValueError: # basically there are not enough fields to split because the "finished"  msg was sent to the sink
		break
		# exit loop-- last read received
	
	mapped_mapperID[counter] = int(mapperID)
	mapped_chrs[counter] = int(chr[3:])
	mapped_pos[counter] = int(position)
	mapped_perfect[counter] = edit == "0"
	
	total_mappings +=1
	counter +=1
	
	if counter == MAX_READS:
		requester.send("SINK UPDATE " + str(counter))
		_ = requester.recv()
		
		#write to hdf5 file and reset arrays
		out_data = np.empty([len(mapped_pos)],dtype='u4,b')
		out_data['f0'] = mapped_pos
		out_data['f1'] = mapped_perfect
		
		# store into temporary HDF5 file.
		for chr in range(1,25):
			mask = mapped_chrs==chr
			tables[chr-1].append(out_data[mask])
		
		for chr in range(1,25):
			mask = mapped_chrs==chr+30
			tables[24+chr-1].append(out_data[mask])
		
		#requester.send("SINK UPDATE DONE")
		#_ = requester.recv()
		
		mapped_chrs = np.zeros([MAX_READS],dtype=np.uint8)
		mapped_pos = np.zeros([MAX_READS],dtype=np.uint32)
		mapped_perfect = np.zeros([MAX_READS],dtype=np.bool)
		mapped_mapperID = np.zeros([MAX_READS],dtype=np.uint8)
		
		counter = 0
		h5file.flush()


# FINAL WRITE OUT

requester.send("SINK UPDATE " + str(counter))
_ = requester.recv()

#write to hdf5 file and reset arrays
out_data = np.empty([len(mapped_pos)],dtype='u4,b')
out_data['f0'] = mapped_pos
out_data['f1'] = mapped_perfect

# store into temporary HDF5 file.
for chr in range(1,25):
	mask = mapped_chrs==chr
	tables[chr-1].append(out_data[mask])

for chr in range(1,25):
	mask = mapped_chrs==chr+30
	tables[24+chr-1].append(out_data[mask])


h5file.flush()
h5file.close()

requester.send("SINK FINISHED " + str(total_mappings))
_ = requester.recv()