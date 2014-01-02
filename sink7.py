# required modules: 
# module load zlib hdf5/1.8.4 
# module load hdf5/1.8.4 (or 1.8.2 on some nodes)

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

print sys.version

import zmq

def errorMsg(msg):
	requester.send(msg)
	_ = requester.recv()
	exit(0)

context = zmq.Context(8)


controlleraddress = sys.argv[1]
controllerport = int(sys.argv[2])

print "Connecting to controller @ " +str(controlleraddress) + ":" + str(controllerport)

MAPPER2SINK_PORT = controllerport +3 #formerly  :6000

print "Will be using port %d as mapper-->sink port" % MAPPER2SINK_PORT



requester = context.socket(zmq.REQ);
requester.connect("tcp://" + controlleraddress + ":" + str(controllerport))

nodeaddress = os.uname()[1]

requester.send("SINK REGISTER " + nodeaddress)
_ = requester.recv()

try:
	insock = context.socket(zmq.PULL)
	insock.bind("tcp://*:"+ str(MAPPER2SINK_PORT))
except zmq.core.error.ZMQError:
	# port already in use!
	errorMsg("ERROR SINK mapper insock port not available!")

quit_signal = False

total_mappings = 0

requester.send("SINK GETSAMPLEINFO")
msg = requester.recv()
hdf5_outfile, sampleID = msg.split(" ")
print "outfilename: " , hdf5_outfile
print "sampleID: " , sampleID



requester.send("SINK GETTRANSLATETABLE")
translate_table_fn = requester.recv()
print "translate table received: " , translate_table_fn



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
temp_hdf5_filename = "/var/tmp/"+str(sampleID)+".tmp.h5"

try:
	h5file = openFile(temp_hdf5_filename, mode = "w", title = "Temp")
except IOError:
	requester.send("ERROR SINK cannot open temp HDF5 file: " + temp_hdf5_filename)
	_ = requester.recv()
	exit(0)

group = h5file.createGroup("/",'mappings','Mappings')


try:
	h5file_out = openFile(hdf5_outfile, mode = "w", title = "Sample " + str(sampleID))
except IOError:
	requester.send("ERROR SINK cannot open final HDF5 file: " + hdf5_outfile)
	_ = requester.recv()
	exit(0)


class mapped_read(IsDescription):
	pos = UInt32Col(pos=0)
	perfect_match = BoolCol(pos=1)

tables = []
for chr in range(1,25):
	tables.append(h5file.createTable(group,"chr"+str(chr),mapped_read,"chr"+str(chr)))#,expectedrows=mapping_exprows[chr]))

h5file.flush()

while 1: 
	data_msg = insock.recv()
	try:
		mapperID, chr, position, edit = data_msg.split(" ")
	except ValueError: # basically there are not enough fields to split because the "finished"  msg was sent to the sink
		# a finish command was received from the mappers
		
		mapperID, mappingCnt , mappedSeqCnt = data_msg.split(" ")
		print "SINK MAPPERFINISH %s %s %s" % (mapperID, mappingCnt , mappedSeqCnt)
		sys.stdout.flush()
		
		requester.send("SINK MAPPERFINISH %s %s %s" % (mapperID, mappingCnt , mappedSeqCnt))
		cmd = requester.recv()
		
		# send received # of reads to controller		
		if cmd == 'SINK KILL':
			break
			# exit loop-- all mappers have finished now and the last read was received
		else:
			continue
		
	
	#mapped_mapperID[counter] = int(mapperID)
	mapped_chrs[counter] = int(chr[3:])
	mapped_pos[counter] = int(position)
	mapped_perfect[counter] = edit == "0"
	
	
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
		
		
		#requester.send("SINK UPDATE DONE")
		#_ = requester.recv()
		
		mapped_chrs = np.zeros([MAX_READS],dtype=np.uint8)
		mapped_pos = np.zeros([MAX_READS],dtype=np.uint32)
		mapped_perfect = np.zeros([MAX_READS],dtype=np.bool)
		mapped_mapperID = np.zeros([MAX_READS],dtype=np.uint8)
		total_mappings += counter
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


h5file.flush()

##### receive contig read count data!!!

total_mappings += counter

requester.send("SINK FINISHED " + str(total_mappings))
contig_data = requester.recv()
print contig_data


##################  CLOSE PORTS

requester.close()
insock.close()

context.term()

################ POST PROCESSING ######################

print "Grouping Reads and Translating from hg19_exome to hg19"
t_start = time.time()
#translate_table_ccds = '/net/grc/shared/scratch/nkrumm/exome_CCDS.txt'
#translate_table_refseq = '/net/grc/shared/scratch/nkrumm/exome_ref.txt'
#translate_table_fn = '/net/grc/shared/scratch/nkrumm/BROAD/broad_ESP.targets.regions.merged.offset_table.masked.txt'


# HDF5 row descriptors
class position(IsDescription):
    pos = UInt32Col(pos=0)
    read_starts = UInt16Col(pos=1)
    read_starts_pm = UInt16Col(pos=2)

class read_info(IsDescription):
	contig = StringCol(10,pos=0)
	unique_reads = UInt32Col()
	mapped_locations = UInt32Col()



def readTranslateTable(tt_file_name):
	tt_starts = [{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}] #list of chromosome dictionaries for the starts 
	tt_stops =  [{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}] #and for the stops of each exon.
	f = open(tt_file_name)
	for line in f: #read in file.
		line = line.strip("\n").split("\t")
		chr = int(line[0][3:])
		tt_starts[chr-1][int(line[4])] = int(line[1])
		tt_stops[chr-1][int(line[5])] = int(line[2])
	return tt_starts,tt_stops

def groupAndCountReads(values,mask=None):
	if mask is None:
		values = np.sort(values) #sort reads
		diff = np.concatenate(([1],np.diff(values)))
		idx = np.concatenate((np.where(diff)[0],[len(values)]))
		return values[idx[:-1]],np.diff(idx) # read positions, counts
	else:
		values = values[mask]
		values = np.sort(values)
		diff = np.concatenate(([1],np.diff(values)))
		idx = np.concatenate((np.where(diff)[0],[len(values)]))
		return values[idx[:-1]],np.diff(idx) # read positions, counts

def translateCoords(a,chr,tt_starts,tt_ix):
	k = np.searchsorted(tt_ix[chr-1],a) - 1 # find the index in the translate indices generated above
	offset = [tt_starts[chr-1][kk] - kk for kk in tt_ix[chr-1][k]]
	return a+offset

tt_starts,tt_stops = readTranslateTable(translate_table_fn)
#tt_starts_refseq,tt_stops_refseq = readTranslateTable(translate_table_refseq)

#make index arrays. These are simply for numpy speedup in translating speed.
ix_exons = {}
for i in range(0,24):
	ix_exons[i] = np.sort(tt_starts[i].keys())


#### Set up outfile for final hdf5 read depth file

group = h5file_out.createGroup("/",'read_starts','Read Starts')
outtables = []
for chr in range(1,25):
	outtables.append(h5file_out.createTable(group,"chr"+str(chr),position,"chr"+str(chr)))#,expectedrows=rs_exprows[chr]))

info_group = h5file_out.createGroup("/",'info_group','Read Info')
info_table = h5file_out.createTable(info_group,"info_table",read_info,"Read Info Table")


#### add contig read count information

contigs = contig_data.split(";")
for c in contigs:
	if len(c) > 0:
		name, readcount = c.split(":")
		m = info_table.row
		m["unique_reads"] = int(readcount)
		m["contig"] = name
		m.append()


### sort and group reads
for chr in range(1,25):
	in_table = h5file.root.mappings._f_getChild("chr" + str(chr))
	
	out_table = h5file_out.root.read_starts._f_getChild("chr" + str(chr))
	
	t1 = time.time()
	values = in_table.read(field='pos')
	perfect_match = in_table.read(field='perfect_match')
	#print chr, len(values)
	
	if len(values) > 0:
		a,out_count = groupAndCountReads(values)
		out_position   = translateCoords(a,chr,tt_starts,ix_exons)
		ix = np.argsort(out_position) # get final sort		
		if np.any(perfect_match):
			a,out_count_pm = groupAndCountReads(values,mask=perfect_match)
			out_position_pm   = translateCoords(a,chr,tt_starts,ix_exons)
			ix_pm = np.argsort(out_position_pm) # get final sort
			out_count_pm_matched = np.zeros(len(out_count))
			out_count_pm_matched[np.searchsorted(out_position[ix],out_position_pm[ix_pm])] = out_count_pm[ix_pm]
		else:
			out_count_pm_matched = np.zeros(len(out_count))
		
		out_data = np.empty([len(out_position)],dtype='u4,u2,u2')
		out_data['f0'] = out_position[ix]
		out_data['f1'] = out_count[ix]
		out_data['f2'] = out_count_pm_matched
	else:
		out_data = np.empty([0],dtype='u4,u2,u2')
		out_data['f0'] = np.array([])
		out_data['f1'] = np.array([])
		out_data['f2'] = np.array([])

	t2 = time.time()
	
	
	# write to file
	
	out_table.append(out_data)
	print "done with chr", chr, "grouping time: ", t2-t1,"writing time: ", time.time()- t2

h5file_out.flush()
h5file_out.close()
h5file.close()

os.system("rm -f " + temp_hdf5_filename)