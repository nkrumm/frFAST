import pysam
import time

sampleID = '3616'
source_filename ='/net/grc/shared/released_exomes/'+sampleID+'/'+sampleID+'.merged.sorted.nodups.realigned.all_reads.bam'

s = pysam.Samfile(source_filename,'rb')


#s.seek(0)
reads = s.fetch(until_eof=True)
#reads = s.fetch()
t1 = time.time()
cnt = 0
while cnt < 1000000:
	r = reads.next()
	seq = reverseComplement(r.seq)
	
	cnt += 1

t2 = time.time()

print t2-t1

########################

from subprocess import Popen, PIPE
import string
import time

def samline(fname): 
	f = Popen(['/Users/nkrumm/Downloads/samtools-0.1.18/samtools', 'view', fname], stdout=PIPE)
	for line in f.stdout:
		yield line

def samlines(fname): 
	f = Popen(['samtools', 'view', fname], stdout=PIPE)
	while 1:
		a = f.stdout.readlines(1000000) # 1000000 = 2400 lines, 10000000 = 24000
		for i in range(len(a)):
			yield a[i]

complement = string.maketrans('atcgnATCGN', 'tagcnTAGCN')

def reverseComplement(sequence):
	return sequence.translate(complement)[::-1]

source_filename = '/net/eichler/vol8/home/nkrumm/EXOMES/3615.merged.sorted.nodups.realigned.all_reads.bam'
samfile = samline(source_filename)

rc_flag = 16
number_splits =1 
split_length = 36
read_length = 50
rc_offset_start = read_length - split_length
read_counter = 0
total_reads = 0
t1 = time.time()
cnt = 0
while read_counter < 1000000:
	try:
		read = samfile.next().split("\t",10)
	except StopIteration:
		endOfFile = True
		break
 	
 	if read[1] == '83' or read[1] == '147':# needs to be reverse complemented.
 		msg = "%s %s" % (read[0], reverseComplement(read[9][rc_offset_start:]))
 	else:
		msg = "%s %s" % (read[0],read[9][0:split_length])
	
	read_counter += 1

t2 = time.time()

print t2-t1
