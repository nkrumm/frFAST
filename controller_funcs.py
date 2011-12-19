import zmq
import time
import operator
import numpy as np
import tempfile
import os
import subprocess
import sys
import signal
import curses
import mapper

def scanPorts():
	import subprocess
	import re	
	cmd = "/opt/c3-4/cexec netstat -vatn"
	args =  cmd.split(" ")
	output,error = subprocess.Popen(args,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
	ports = []
	for line in output.split("\n"):
		if line[0:3] == 'tcp':
			for rexp in re.findall(":[0-9]+", line):
				ports.append(int(rexp.strip(":")))
	
	return set(ports)

def writeQSubFile(command):
	tf = tempfile.NamedTemporaryFile('w')
	tf.write("""
	source ~/.bash_profile
	module load modules modules-init modules-gs modules-eichler
	module load python/2.7.2
	""" + command)
	tf.flush()
	return tf

def submitQSub(command):
	args=command.split(" ")
	output,error = subprocess.Popen(args,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
	return output, error


def deleteQSub(jobID):
	command = "qdel " + str(jobID)
	args=command.split(" ")
	output,error = subprocess.Popen(args,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
	return output, error

def logMessage(f_log, logLevel, msg):
	if logLevel != 0:
		f_log.write(msg + "\n")

def updateScreen(screen,mappers):
	if screen != None:
		screen.clear()
		screen.box()
		screen.addstr(0,3,"   MAPPING JOBS   ")
		screen.addstr(2,3,"mapperID\tnode\tstatus\ttask\trank")
		row = 0
		for _mapper in mappers.values():
			wstr = "Working" if _mapper["working"] else "Free"
			screen.addstr(3+row,3,str(_mapper["mapperID"]) + "\t\t" + str(_mapper["node"]) + "\t" + wstr + "\t" + str(_mapper["task"]) + "\t" + str(_mapper["rank"]))
			row += 1
		screen.refresh()

def updateMessages(msgHistory, screen, f_log, logLevel, msg):
	if len(msgHistory) > 15:
		_ = msgHistory.pop(0)
	t = time.localtime()
	timestr = "[%02d:%02d:%02d] " % (t.tm_hour, t.tm_min, t.tm_sec)
	msgHistory.append(timestr + msg)
	logMessage(f_log, logLevel, timestr + msg)
	if screen != None:
		screen.clear()
		screen.box()
		screen.addstr(0,3,"   MESSAGES   ")
		row = 0
		for m in msgHistory:
			screen.addstr(2+row,3,m)
			row += 1
		screen.refresh()


def updateStats(statsData, start_time, screen, data): # data is a dictionary! {"var" : value}
	# update global data storage
	for d in data.keys():
		statsData[d]["value"] += int(data[d])
	if screen != None:
		screen.clear()
		screen.box()
		screen.addstr(0,3,"   STATISTICS   ")
		row = 0
		
		for d in statsData.keys():
			if d != 'contigreadcnt':
				screen.addstr(2+row,3,statsData[d]["desc"] + ":\t" + str(statsData[d]["value"]))
				row += 1
		
		# update time
		elapsed_sec = time.time() - start_time 
		if elapsed_sec > 60:
			screen.addstr(3+row,3,"Elapsed Time: "  + str(int(elapsed_sec/60)) + "m:" + str(int(elapsed_sec % 60)) + "s")
		else:
			screen.addstr(3+row,3,"Elapsed Time: 0m:"  + str(int(elapsed_sec)) + "s")
		screen.refresh()

def quitController(quitMessage, logDirectory, job_ids, f_log, f_summary, f_contigs, gui):
	# an error occurred-- mark and start on the next one
	f_error = open(logDirectory + "/error.log",'w')
	if quitmessage != None:
		f_error.write("ERROR: timeout exceeded!\n")
	f_error.write('EXITING! Killing qsub jobs:\n')
	for jobID in job_ids.values():
		o,e = deleteQSub(jobID)
		f_error.write(o)
	
	f_error.close()
	f_log.close()
	f_summary.close()
	f_contigs.close()
	
	if gui:
		curses.endwin()
	sys.exit(0)
