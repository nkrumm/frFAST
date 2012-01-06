
file = '/net/eichler/vol16/a5ko/Exomes/sequences/11013.mo/sequences/11013.mo/lib1/fastq_split/11013.mo.5286.merged.sorted.nodups.realigned.all_reads.bam.1_part0.fastq.gz'


import zlib
import time
import copy
import pysam
def zlib_lines(f, d, prev_fieldcnt = 0, prev_tailstring = ''):
	# http://stackoverflow.com/questions/2423866/python-decompressing-gzip-chunk-by-chunk
	CHUNKSIZE=1024*32
	buffer=f.read(CHUNKSIZE)
	read = pysam.AlignedRead()
	#fieldcnt = 0
	#prev_tailstring = ''
	while buffer:
		outstr = d.decompress(buffer)
		split_lines = outstr.split("\n")
		
		if prev_tailstring != '':
			split_lines[0] = prev_tailstring+split_lines[0]
			#print "starting string: %s " % split_lines[0]
		#else:
			#print "no change to starting string: %s " % split_lines[0]
		
		if outstr[-1] != "\n":
			i = iter(split_lines[:-1])
			tailstring = split_lines[-1]
			#print "tailstring: %s " % tailstring
		else:
			i = iter(split_lines[:-1])
			tailstring = ''
			#print "tailstring: %s " % tailstring
			#print "outstr[-1] = : %s" % outstr.split("\n")[-1]
		
		if prev_fieldcnt not in [0,4]: #the previous chunk ended in the middle of a read
			#print "fieldcount is: %d" % fieldcnt
			if prev_fieldcnt == 1:
				read.seq = i.next()
				_ = i.next()
				_ = i.next()
			elif prev_fieldcnt == 2:
				_ = i.next()
				_ = i.next()
			elif prev_fieldcnt == 3:
				_ = i.next()
			prev_fieldcnt = 0
			#print read.qname, read.seq
			yield read, prev_fieldcnt, prev_tailstring
			
			
		while 1:
			fieldcnt = 0
			try:
				read.qname = i.next()
				fieldcnt += 1
				
				read.seq = i.next()
				fieldcnt += 1
				
				_ = i.next() # iterate past the '+' and the ...
				fieldcnt += 1
				
				_ = i.next() # ... quality string
				fieldcnt +=1
				
				yield read, prev_fieldcnt, prev_tailstring
				
			except StopIteration:
				# the end of the chunk was hit.
				break
			
		
		# read in the next buffer chunk
		buffer=f.read(CHUNKSIZE)
		prev_tailstring = tailstring
		prev_fieldcnt = fieldcnt


del f
f = open(file,'rb')
c = 0
d = zlib.decompressobj(16+zlib.MAX_WBITS)

source = zlib_lines(f,d)
seek_pos = []
reads = []
dobjs = []
fieldcnts = []
tailstrings = []
t1 = time.time()
while 1:
	
	try:
		if c % 4000 == 0:
			read,fieldcnt,tailstring = source.next()
			seek_pos.append(f.tell())
			dobjs.append(d.copy())
			reads.append(read.seq)
			fieldcnts.append(fieldcnt)
			tailstrings.append(tailstring)
			#print c, fieldcnt
		else:
			_ = source.next()	
		c+=1
	except StopIteration:
		break

t2 = time.time()
print "Done in ", str(t2-t1), ", lines: ", str(c)

f.seek(seek_pos[5])
source = zlib_lines(f,dobjs[5],prev_fieldcnt=fieldcnts[5],prev_tailstring=tailstrings[5])
print source.next()[0]

seq = ''
while seq != reads[10]:
	seq = source.next().seq
	print seq


f.seek(seek_pos[5])
source = zlib_lines(f,dobjs[5])
print source.next()

############################



from subprocess import Popen, PIPE
import gzip
import time


def pythongziplines(fname): 
	f = gzip.open(fname, 'rb')
	while f.readline():
		yield line

source = pythongziplines(file)

c = 0
t1 = time.time()
while 1:
	try:
		_ = source.next()
		c+=1
	except StopIteration:
		break

t2 = time.time()
print "Done in ", str(t2-t1), ", lines: ", str(c)
f.close()





##############################################################

def gziplines(fname): 
	f = Popen(['zcat', fname], stdout=PIPE)
	for line in f.stdout:
		yield line


source = gziplines(file)

c = 0
t1 = time.time()
while 1:
	try:
		_ = source.next()
		c+=1
	except StopIteration:
		break

t2 = time.time()
print "Done in ", str(t2-t1), ", lines: ", str(c)
f.close()




##############################################################
import shutil
import cStringIO as StringIO


c = 0
t1 = time.time()
output = StringIO.StringIO()
shutil.copyfileobj(gzip.open(file,'rb'), output)
output.seek(0)

while 1:
	try:
		_ = output.next()
		c+=1
	except StopIteration:
		break

t2 = time.time()
print "Done in ", str(t2-t1), ", lines: ", str(c)
f.close()




###################################
import subprocess
c = 0
t1 = time.time()
p = subprocess.Popen(["zcat", file], stdout = subprocess.PIPE)
f = cStringIO.StringIO(p.communicate()[0])

seek_pos = []
for line in f:
	if c % 4000 == 0:
		seek_pos.append(f.tell())
	c+=1

t2 = time.time()
print "Done in ", str(t2-t1), ", lines: ", str(c)