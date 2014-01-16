# module load python/2.7.2
import argparse
import zmq
import time
import os
import subprocess
import sys
import signal
import curses
import mapper
from controller_funcs import scanPorts, writeQSubFile, submitQSub, deleteQSub, logMessage, updateScreen, updateMessages, updateStats, quitController

start_time = time.time()
################################################################################################
                               ######### OPTIONS ###########
################################################################################################
parser = argparse.ArgumentParser(description='This is frFAST (frankenFAST), a lightweight pipeline designed to calculate read-depth across an exome or genome using the mrsFAST mapper. \nNik Krumm, 2011')
parser.add_argument('--source', required=True, metavar='/path/to/source_file.bam ', type=str, nargs='?',help="Source to be bam file to be mapped.")
parser.add_argument('--output', required=True, metavar='/path/to/output_file.h5', type=str, nargs='?',help="Output location of HDF5 file")
parser.add_argument('--log_directory','--log_dir', required=True, metavar='/path/to/log_directory', type=str, nargs='?',help="Location of log directory. Will be created if necessary.")
parser.add_argument('--error_log', required=False, default="", metavar='/path/to/universal_error_file', type=str, nargs='?',help="Location of central error file. Output is appended to file, but file is created if necessary.")
parser.add_argument('--sampleID', required=True, metavar='sampleID', type=str, nargs='?',help="Unique sampleID string")
parser.add_argument('--index_dir', metavar='/path/to/index/folder', type=str, default="/net/grc/shared/scratch/nkrumm/INDEX/default_exome", nargs='?',help="Location of folder containing index file and index.fa which will be copied to each mapping node via rsync.")
parser.add_argument('--dont-rsync-index', action="store_true", default=False, help="Do not use rsync to copy mrsFAST index to node (--index_dir should point to accessible index in this case)")
parser.add_argument('--index', metavar='index.fa', type=str, default="default_exome.fa", nargs='?',help="Name of the index.fa file used for mapping, must be in folder specified by index_dir \n Default: default_exome.fa")
parser.add_argument('--translate_table', metavar='/path/to/translate_table.txt', type=str, default="/net/grc/shared/scratch/nkrumm/translate_tables/default_exome.translate.txt", nargs='?',help="Path to correct translate table used to remap exome-based coordinates to genome based coordinates. Default is based on Nimblegen V2 capture.")
parser.add_argument('--port','-p', metavar='8000', type=int, nargs="?", default = 8000,\
    help="TCP port offset to use. Will use FOUR CONSECUTIVE ports starting with this value. Default port range: 8000-8004")
parser.add_argument('--chunksize', metavar='1000000', type=int, default=1000000, help='Number of reads to map per chunk.')
parser.add_argument("--disable-port-scan", '--disable_port_scan',action='store_true',\
    help="Disable initial port scan. Use with caution and only with ports you know are open!")
parser.add_argument("--disable-gui", '--disable_gui',action='store_true',\
    help="Disable dashboard GUI for running in batch/SGE")
parser.add_argument('--timeout','-t', metavar='3600', type=int, nargs="?", default = 0,\
    help="Maximum idle timeout before aborting current sample. Default is no timeout. SGE queue time is not included in this.")
parser.add_argument('--single-host',action='store_true',\
    help="Run using localhost/tcp-loopback only; Disables SGE queue submission and runs mappers/sinks/vents using subprocesses!")
parser.add_argument('--n-mappers',default=8, help="Number of mapping jobs") 

args = parser.parse_args()

if args.disable_gui:
    gui = False
else:
    gui = True
    
sampleID = args.sampleID #str(sys.argv[1])
tcpPortIndex = args.port #int(sys.argv[2])
source_filename = args.source #'/net/grc/shared/released_exomes/'+sampleID+'/'+sampleID+'.merged.sorted.nodups.realigned.all_reads.bam'
sink_outfile = args.output #"/net/grc/shared/scratch/nkrumm/ESP2000/HDF5/" + sampleID + ".h5"

split_length = 36 # only 36 is currently supported!
read_length = 'AUTO'

# set up index location information
if args.dont_rsync_index:
    indexFile = args.index
    index_on_node_location = os.path.join(args.index_dir, args.index)
    INDEXDIRPATH = ""
else:
    indexFile = args.index #"/var/tmp/exome/exome_19_masked.fa"
    INDEXDIRPATH = "/" + args.index_dir.strip("/")
    folder = INDEXDIRPATH.split("/")[-1]
    index_on_node_location = "/var/tmp/" + folder + "/" + indexFile

# translation table to be sent to sink
translate_table_fn = args.translate_table

logLevel = 1 # 1: log messages to file,  0: no logging TODO: turn off/on SGE logs.
logDirectory = args.log_directory + "/" #"/net/grc/shared/scratch/nkrumm/ESP2000/Logs/" + sampleID + "/"

job_names = {"map": "map" + sampleID, "vent": "vent" + sampleID, "sink": "sink" + sampleID}
mem_requested = {"mapper": "1500M", "vent": "2G", "sink": "3G"}

BASEFILE =  os.path.realpath(__file__)
BASEPATH =  os.path.dirname(BASEFILE)

scriptFiles= {"map":  BASEPATH + "/mrsfast_wrapper.py", "vent": BASEPATH + "/ventilator_failsafe3.py", "sink": BASEPATH + "/sink7.py"}
#scriptFiles= {"map":  BASEPATH + "/mrsfast_wrapper.py", "vent": BASEPATH + "/ventilator_downsampler.py", "sink": BASEPATH + "/sink7.py"}
totalMappers = args.n_mappers
maximum_reads_per_destination = args.chunksize


################################################################################################
################################################################################################


# Set up output log files
if os.path.exists(logDirectory) != True:
    os.mkdir(logDirectory)

if args.error_log !="":
    if os.path.exists(args.error_log):
        f_univerror = open(args.error_log,'a')
    else:
        f_univerror = open(args.error_log,'w')
else:
    f_univerror = None

f_log = open(logDirectory + "/messages.log",'w')
f_summary = open(logDirectory + "/summary.log",'w')
f_contigs = open(logDirectory + "/contigs.log",'w')


msgHistory = []
    
statsData = {}
statsData["ventreadcnt"] = {"desc": "Total Ventilated Reads", "value": 0}
statsData["mappedreadcnt"] = {"desc": "Total Mapped Reads", "value": 0}
statsData["sinkmappingcnt"] = {"desc": "Total Mappings Received by Sink", "value": 0}
statsData["totalmappings"] = {"desc": "Total Mappings", "value": 0}
statsData["contigreadcnt"] = {}

############################

if (not args.disable_port_scan) and (not args.single_host):
    print "Scanning ports on cluster..."
    ports = scanPorts()
    
    for i in range(tcpPortIndex, tcpPortIndex+4):
        if i in ports:
            print "Error: ports selected already in use. Please try a different port offset!"
            sys.exit(0)

REQ_REP_PORT = tcpPortIndex # formerly :5555 aka controllerport
PUB_SUB_PORT = tcpPortIndex + 1 # :7000
VENT2MAPPER_PORT = tcpPortIndex + 2 #:8000
MAPPER2SINK_PORT = tcpPortIndex + 3 #:6000

#############################
def SIGINT_handler(signal, frame):
    if args.single_host:
        print 'EXITING with return code 1!'
    else:
        print 'EXITING with return code 1! Killing qsub jobs:'
        for jobID in job_ids.values():
            o,e = deleteQSub(jobID)
            print o
    if gui:
        curses.endwin()
    sys.exit(1)
    #curses.nocbreak(); stdscr.keypad(0); curses.echo()
    #curses.endwin()


# set up controller 
context = zmq.Context()
socket = context.socket(zmq.REP)
socket.bind("tcp://*:" + str(REQ_REP_PORT))

pubsocket = context.socket(zmq.PUB)
pubsocket.bind("tcp://*:" + str(PUB_SUB_PORT))

signal.signal(signal.SIGINT, SIGINT_handler) # handle Ctrl-C gracefully

if args.single_host:
    controllerNodeName = "localhost"
else:
    controllerNodeName = os.uname()[1].split(".")[0] #own hostname
mappers_started = False

job_ids = {}
qsub_cmd = {}


#############################
# START VENT
#
vent_is_ready = False
vent_is_finished = False

if args.single_host:
    qsub_cmd["vent"] = "python %s %s %d" % (scriptFiles["vent"], controllerNodeName, REQ_REP_PORT)
    p = subprocess.Popen(qsub_cmd["vent"], shell=True)
    job_ids["vent"] = p.pid
else:
    tempscript = writeQSubFile("module load samtools/latest python/2.7.2 numpy/1.6.1\n python "+scriptFiles["vent"]+" " + controllerNodeName+" " + str(REQ_REP_PORT))
    qsub_cmd["vent"] = "qsub -pe serial 2 -l h_vmem=%s -l rhel=6 -S /bin/bash -cwd -o %s -e %s -N %s -j y %s" % (mem_requested["vent"], logDirectory, logDirectory,job_names["vent"], tempscript.name)
    output,e = submitQSub(qsub_cmd["vent"])
    job_ids["vent"] = output.split(" ")[2]
    vent_time = 0

#############################
# START SINK
#
sink_is_ready = False
if args.single_host:
    qsub_cmd["sink"] = "python %s %s %d" % (scriptFiles["sink"], controllerNodeName, REQ_REP_PORT)
    p = subprocess.Popen(qsub_cmd["sink"], shell=True)
    job_ids["sink"] = p.pid
else:
    tempscript = writeQSubFile("module load zlib/latest hdf5/1.8.4 python/2.7.2 numpy/1.6.1\n python "+scriptFiles["sink"]+" " + controllerNodeName+" " + str(REQ_REP_PORT))
    qsub_cmd["sink"] = "qsub -pe serial 2 -l h_vmem=%s -l rhel=6 -S /bin/bash -cwd -o %s -e %s -N %s -j y %s" % (mem_requested["sink"], logDirectory, logDirectory,job_names["sink"], tempscript.name)
    output,e = submitQSub(qsub_cmd["sink"])
    job_ids["sink"] = output.split(" ")[2]

####
# SET UP MAPPERS and job details

mappers = mapper.Mappers(totalMappers)
curr_mapperID = 0

ventilated_reads = {}
mapper_received_reads = {}


retry_chunk_queue = []
chunkID_lookup = {}

#############################
# MAIN CONTROL LOOP
#
t_start = time.time()

poller = zmq.Poller()
poller.register(socket, zmq.POLLIN)
poller.register(pubsocket, zmq.POLLOUT)

if gui:
    window = curses.initscr()
    window.border(0)
    window.clear()
    window.refresh()
    screen = curses.newwin(20, 80, 0, 0)
    
    msgScreen = curses.newwin(20, 200, 20, 0)
    statScreen = curses.newwin(10, 120, 10, 80) #h,w,begin_y,begin_x
    infoScreen = curses.newwin(10, 120, 0, 80) #h,w,begin_y,begin_x
    infoScreen.clear()
    infoScreen.box()
    infoScreen.addstr(0,3,"   JOB INFO   ")
    infoScreen.addstr(2,3,"sampleID: " + str(sampleID))
    infoScreen.addstr(3,3,"Input: " + str(source_filename))
    infoScreen.addstr(4,3,"Output: " + str(sink_outfile))
    infoScreen.refresh()
else:
    screen = None
    msgScreen = None
    statScreen = None

updateScreen(screen, mappers.mappers)
updateMessages(msgHistory, msgScreen, f_log, logLevel, "[INFO] Controller running on: " + controllerNodeName)
updateMessages(msgHistory, msgScreen, f_log, logLevel, "[INFO] ZMQ library version: " + str(zmq.__version__))
updateMessages(msgHistory, msgScreen, f_log, logLevel, "[INFO] Will be using index in: " + INDEXDIRPATH)
updateMessages(msgHistory, msgScreen, f_log, logLevel, "[INFO] Using index: " + indexFile)
updateMessages(msgHistory, msgScreen, f_log, logLevel, "[INFO] Using source: " + source_filename)
updateMessages(msgHistory, msgScreen, f_log, logLevel, "[INFO] Index on node is at: " + index_on_node_location)
updateMessages(msgHistory, msgScreen, f_log, logLevel, "[INFO] Setting as sink output file: " + sink_outfile)
updateMessages(msgHistory, msgScreen, f_log, logLevel, "[INFO] Vent script location: " + scriptFiles["vent"])
updateMessages(msgHistory, msgScreen, f_log, logLevel, "[INFO] Submitted vent job: " + qsub_cmd["vent"])
updateMessages(msgHistory, msgScreen, f_log, logLevel, "[INFO] SGE Job id for vent: " + str(job_ids["vent"]))
updateMessages(msgHistory, msgScreen, f_log, logLevel, "[INFO] Sink script location: " + scriptFiles["sink"])
updateMessages(msgHistory, msgScreen, f_log, logLevel, "[INFO] Submitted sink job: " + qsub_cmd["sink"])
updateMessages(msgHistory, msgScreen, f_log, logLevel, "[INFO] SGE Job id for sink: " + str(job_ids["sink"]))
updateMessages(msgHistory, msgScreen, f_log, logLevel, "[INFO] Attempting to register " + str(totalMappers) + " mapping jobs")
updateMessages(msgHistory, msgScreen, f_log, logLevel, "[INFO] mem_requested per mapper set to " + mem_requested["mapper"] + " mapping jobs")

time_last_msg  = time.time()

while True:
    if args.timeout != 0: # timeout enabled
        if time.time() - time_last_msg  >= args.timeout:
            if sink_is_ready and vent_is_ready and mappers_started:
                error_msg = "ERROR [%s] Timeout Exceeded!" % sampleID
                quitController(error_msg, logDirectory, job_ids, f_log, f_summary, f_univerror, f_contigs, gui, args.single_host)
    
    if sink_is_ready and vent_is_ready:
        if not mappers_started:
            # START ZE MAPPERS!
            if args.single_host:
                qsub_cmd["mappers"] = "python %s %s %d %s" % (scriptFiles["map"], controllerNodeName, REQ_REP_PORT, BASEPATH)
                pids = []
                for i in range(mappers.numMappers):
                    p = subprocess.Popen(qsub_cmd["mappers"], shell=True)
                    pids.append(p.pid)
                job_ids["mappers"] = pids
                updateMessages(msgHistory, msgScreen, f_log, logLevel,  "#### STARTED MAPPERS ####")
                #updateMessages(msgHistory, msgScreen, f_log, logLevel,  "job_ids are: %s" % ", ".join(job_ids["mappers"]))
            else:
                if args.dont_rsync_index:
                    tempscript = writeQSubFile("python " + scriptFiles["map"] + " " + controllerNodeName+" " + str(REQ_REP_PORT) + " " + BASEPATH)
                else:
                    tempscript = writeQSubFile("python " + scriptFiles["map"] + " " + controllerNodeName+" " + str(REQ_REP_PORT) + " " + BASEPATH + " --index-dir " + INDEXDIRPATH)
                
                qsub_cmd = "qsub -l h_vmem=%s -l rhel=6 -t 1-%d -S /bin/bash -cwd -o %s -e %s -N %s -j y %s" % (mem_requested["mapper"], mappers.numMappers, logDirectory, logDirectory,job_names["map"],tempscript.name)
                output,e = submitQSub(qsub_cmd)
                #print output,e
                job_array = output.split(" ")[2]
                job_id = job_array.split(".")[0]
                job_ids["mappers"] = [str(job_id)]
                job_name = output.split(" ")[3][2:-2] # get the job name only from '("jobname")'
                
                updateMessages(msgHistory, msgScreen, f_log, logLevel,  "#### STARTED MAPPERS ####")
                updateMessages(msgHistory, msgScreen, f_log, logLevel,  "job_array is: " + str(job_array))
                updateMessages(msgHistory, msgScreen, f_log, logLevel,  "job_id is: " + str(job_id))
                updateMessages(msgHistory, msgScreen, f_log, logLevel,  "job_name is: " + str(job_name))
            
            mappers_started = True
    
    
    socks = dict(poller.poll())
    if socks.get(pubsocket) == zmq.POLLOUT:
        pass
    if socks.get(socket) == zmq.POLLIN:
        message = socket.recv()
        time_last_msg = time.time()
        
        #print message
        try: 
            sender, cmd, data = message.split(" ",2)
        except ValueError: 
            sender, cmd = message.split(" ")
            data = None
        
        #############################
        # VENTILATOR 
        #
        if sender == 'VENT': #process ventilator requests
            if cmd == 'REGISTER':
                ventilator_address = data
                updateMessages(msgHistory, msgScreen, f_log, logLevel, "[VENT] Ventilator registered @ "+ ventilator_address)
                if gui:
                    infoScreen.addstr(5,3,"Ventilator Node: " + str(ventilator_address))
                    infoScreen.refresh()
                socket.send("ok")
            
            elif cmd == 'GETSOURCE':
                updateMessages(msgHistory, msgScreen, f_log, logLevel, "[VENT] Sending read source file to ventilator ("+ source_filename +")")
                socket.send(source_filename)
            
            elif cmd == 'READY':
                updateMessages(msgHistory, msgScreen, f_log, logLevel,  "[VENT] Ventilator ready!")
                socket.send("ok")
                vent_is_ready = True
            
            elif cmd == 'GETDEST':
                
                dest,nextMapperID = mappers.getNextDestination()
                
                if nextMapperID != None:
                    if len(retry_chunk_queue) == 0:
                        chunkID = int(nextMapperID)
                        chunkID_lookup[int(nextMapperID)] = chunkID
                    else:
                        chunkID = int(retry_chunk_queue.pop())
                        chunkID_lookup[int(nextMapperID)] = chunkID
                    
                    msg = "%s %d %d" % (str(dest), nextMapperID, chunkID)
                    socket.send(msg)
                    updateMessages(msgHistory, msgScreen, f_log, logLevel, "[VENT] Sending chunkID = %d to mapperID = %d at node = %s" % (chunkID, nextMapperID, str(dest)))
                    updateMessages(msgHistory, msgScreen, f_log, logLevel, "[VENT] Sending GO to "+ str(nextMapperID))
                    pubsocket.send(str(nextMapperID) + " GO") # TODO: send the port to be used by mapper
                    vent_time = time.time()
                    mappers.updateMapperWorkingStatus(nextMapperID, working=True)
                    mappers.updateMapperTask(nextMapperID, "Ventilating")
                    updateScreen(screen,mappers.mappers)
                else:
                    socket.send("0")
                    #updateMessages(msgHistory, msgScreen, f_log, logLevel, "VENT PING")
            
            elif cmd == 'GETNUMREADS':
                socket.send(str(maximum_reads_per_destination))
            elif cmd == 'GETSPLITPARAMS':
                socket.send(str(split_length) + " " + str(read_length))
                
            elif cmd == 'DONE':
                read_counter, total_reads, destination, mapperID = data.split(" ")
                updateMessages(msgHistory, msgScreen, f_log, logLevel, "[VENT] Ventilator done with writing "+str(read_counter)+ " reads to " + str(mapperID) + "; total time: " + str(time.time()-vent_time))
                
                mappers.updateMapperTask(mapperID, "Waiting")
                #mappers.printMappers()
                updateScreen(screen,mappers.mappers)
                updateStats(statsData, start_time, statScreen,{"ventreadcnt":read_counter})
                ventilated_reads[mapperID] = int(read_counter)
                
                if mapper_received_reads.has_key(mapperID) and ventilated_reads[mapperID] != mapper_received_reads[mapperID]:
                    updateMessages(msgHistory, msgScreen, f_log, logLevel, "[MAPPER] ERROR! mapper did not receive all the ventilated reads! " + data)
                    
                    # the mapper is waiting here. mapper needs to restart, so tell it via subsocket
                    pubsocket.send(str(mapperID) + " RESTART")
                    
                    socket.send("RESTART") # this goes back to ventilator! NOT to the mapper
                    mappers.stopMapper(int(mapperID))
                    # add chunkID to retry queue for ventilator
                    retry_chunk_queue.append(chunkID_lookup[int(mapperID)])
                    
                else: #everything is OK, the mapper should proceed
                    socket.send("ok")
                    pubsocket.send(str(mapperID) + " GO")
            
            elif cmd == 'FINISHED':
                read_counter, total_reads = data.split(" ")
                updateMessages(msgHistory, msgScreen, f_log, logLevel, "[VENT] Ventilator done with all reads!; total reads ventilated = " + str(total_reads))
                socket.send("ok")
                vent_is_finished = True
                if mappers.getNumWorking() == 0: #this is the exception case where the mappers finished before the vent
                    # send a dummy MAPPERFINISH signal to the sink via its pull socket
                    updateMessages(msgHistory, msgScreen, f_log, logLevel, "SENDING DUMMY SIGNAL")
                    tempsocket = context.socket(zmq.PUSH)
                    tempsocket.connect("tcp://" + sinkaddress + ":" + str(MAPPER2SINK_PORT))
                    tempsocket.send("0 0 0")
                    tempsocket.close()
                    updateMessages(msgHistory, msgScreen, f_log, logLevel, "... SUCCESSFUL")
                    
                    # return the cmd with SINK KILL, as in the other section
                    tempmsg = socket.recv()
                    updateMessages(msgHistory, msgScreen, f_log, logLevel, "got back: " + str(tempmsg))
                    socket.send("SINK KILL")
            else:
                pass
        
        ##################
        # MAPPER
        #   
        elif sender == 'MAPPER':
            if cmd == 'REGISTER':
                # should get the hostname here and register
                mapperID = mappers.registerMapper(data) # registerMapper(data) 
                socket.send(str(mapperID))
            elif cmd == 'GETINDEX':
                socket.send(index_on_node_location)
            elif cmd == 'GETSINKADDRESS':
                socket.send("tcp://" + sinkaddress + ":" + str(MAPPER2SINK_PORT))
            elif cmd == 'GETVENTADDRESS':
                socket.send("tcp://" + ventilator_address + ":" + str(VENT2MAPPER_PORT))
            elif cmd == 'INPUT':
                mapperID, inputCnt, discardedCnt = data.split(" ")
                if int(inputCnt) == 0:
                    # All the reads were discarded!
                    # the mapper will restart, but we need to delete the mappers entry here, as this mapperID is now dead!
                    updateMessages(msgHistory, msgScreen, f_log, logLevel, "[MAPPER] No valid sequences in received reads! " + data)
                    mappers.updateMapperTask(int(mapperID), "Done")
                    mappers.stopMapper(int(mapperID))
                    socket.send("RESTART")
                    
                else:
                    mappers.updateMapperTask(mapperID, "Got Input")
                    mapper_received_reads[mapperID] = int(inputCnt) + int(discardedCnt)
                    if ventilated_reads.has_key(mapperID) and ventilated_reads[mapperID] != mapper_received_reads[mapperID]:
                        updateMessages(msgHistory, msgScreen, f_log, logLevel, "[MAPPER] ERROR! mapper did not receive all the ventilated reads! " + data)
                        # mapper needs to restart
                        socket.send("RESTART")
                        mappers.stopMapper(int(mapperID))
                        # add chunkID to retry queue for ventilator
                        retry_chunk_queue.append(chunkID_lookup[int(mapperID)])
                        
                    elif ventilated_reads.has_key(mapperID) == False: #the mapper should WAIT until ventilator is 
                        updateMessages(msgHistory, msgScreen, f_log, logLevel, "[MAPPER] %s Waiting" % mapperID)
                        socket.send("WAIT")
                    else: # everything is OK and the mapper can proceed
                        socket.send("OK")
            
            elif cmd == 'UPDATE':
                #updateMessages(msgHistory, msgScreen, f_log, logLevel, "UPDATE received: "+ data)
                mapperID, contig, mappingCnt, mappedSeqCnt = data.split(" ")
                mappers.updateMapper(int(mapperID), contig, int(mappingCnt), int(mappedSeqCnt))
                mappers.updateMapperTask(mapperID, "Mapping " + str(contig))
                totalReads = sum(statsData["contigreadcnt"][mapperID].values())
                statsData["contigreadcnt"][mapperID][contig] = int(mappedSeqCnt) - totalReads
                updateScreen(screen,mappers.mappers)
                socket.send("ok")
            elif cmd == 'START':
                updateMessages(msgHistory, msgScreen, f_log, logLevel, "[MAPPER] START received: "+ data)
                mapperID, seqListSize, mappedSeqCnt = data.split(" ")
                mappers.startMapper(int(mapperID), int(seqListSize), int(mappedSeqCnt))
                mappers.updateMapperTask(mapperID, "Mapping")
                socket.send("ok")
                statsData["contigreadcnt"][mapperID] = {}
            elif cmd == 'FINISH':
                updateMessages(msgHistory, msgScreen, f_log, logLevel, "[MAPPER] FINISH received: "+ data)
                socket.send("ok")
                updateStats(statsData, start_time, statScreen,{"mappedreadcnt":mappedSeqCnt, "totalmappings":mappingCnt})
                
            else:
                pass
        
        elif sender == 'SINK':
            if cmd == 'UPDATE':
                t2 = time.time()
                #updateMessages(msgHistory, msgScreen, f_log, logLevel, "[" + str(t2-t_start) + "] " + data)
                updateStats(statsData, start_time, statScreen,{"sinkmappingcnt":int(data)})
                socket.send("ok")
            elif cmd == 'REGISTER':
                sinkaddress = data
                updateMessages(msgHistory, msgScreen, f_log, logLevel, "[SINK] Sink registered @ "+ sinkaddress)
                socket.send("ok")
                if gui:
                    infoScreen.addstr(6,3,"Sink Node: " + str(sinkaddress))
                    infoScreen.refresh()                
            elif cmd == 'GETSAMPLEINFO':
                updateMessages(msgHistory, msgScreen, f_log, logLevel, "[SINK] Sending sample info to sink: " + sink_outfile + " " + sampleID)
                socket.send(sink_outfile + " " + sampleID)
            elif cmd == 'GETTRANSLATETABLE':
                updateMessages(msgHistory, msgScreen, f_log, logLevel, "[SINK] Sending translate table to sink: " + translate_table_fn)
                socket.send(translate_table_fn)
            elif cmd == 'READY':
                updateMessages(msgHistory, msgScreen, f_log, logLevel, "[SINK] Sink ready!")
                socket.send("ok")
                sink_is_ready = True
            elif cmd == 'MAPPERFINISH':
                if int(mapperID) != 0: #zero is a special exception, received only when SINK KILL signal is imminent.
                    updateMessages(msgHistory, msgScreen, f_log, logLevel, "[SINK] MAPPERFINISH received: "+ data)
                    mapperID, mappingCnt, mappedSeqCnt  = data.split(" ")
                    mappers.updateMapperTask(int(mapperID), "Done")
                    mappers.stopMapper(int(mapperID), int(mappingCnt), int(mappedSeqCnt))
                    
                if mappers.getNumWorking() == 0:
                    if vent_is_finished:
                        updateMessages(msgHistory, msgScreen, f_log, logLevel, "ALL MAPPERS DONE!")
                        time.sleep(10)
                        socket.send("SINK KILL")
                    else:
                        # the mappers finished before the vent could tell the controller it is done. caveat: if the vent is NOT done but the NumWorking value is 0, this is also triggered
                        # the kill signal will be sent when the ventilator sends its finish signal
                        socket.send("ok")
                else:
                    socket.send("ok")
            
            
            elif cmd == 'FINISHED':
                updateMessages(msgHistory, msgScreen, f_log, logLevel, data)
                #updateStats(statsData, start_time, statScreen,{"sinkmappingcnt":int(data)})
                
                import pickle
                f=open(logDirectory + "/statsData.pickle",'w')
                pickle.dump(statsData,f)
                f.close()
                
                contigs = {}
                # header line
                f_contigs.write('mapperID\t' + '\t'.join(sorted(statsData["contigreadcnt"].values()[0].iterkeys())) + "\n")
                
                for mapperID in sorted(statsData["contigreadcnt"].iterkeys()):
                    f_contigs.write(str(mapperID) + "\t")
                    for contig in sorted(statsData["contigreadcnt"][mapperID].iterkeys()):
                        readcnt = statsData["contigreadcnt"][mapperID][contig]
                        f_contigs.write(str(readcnt) + "\t")
                        try: contigs[contig] += readcnt
                        except KeyError: contigs[contig] = readcnt
                    f_contigs.write("\n")
                
                
                #summary line
                f_contigs.write("totals\t")
                f_contigs.write('\t'.join([str(contigs[a]) for a in sorted(contigs.iterkeys())]))
                
                # send final read count information to sink for final incorporation into HDF5 file
                contig_data = ';'.join([str(a) +":" +str(b) for a,b in zip(contigs.keys(),contigs.values())])
                socket.send(contig_data)
                
                f_contigs.close()
                
                f_summary.write("Source\t" + str(statsData["ventreadcnt"]["desc"]) + "\t" + str(statsData["mappedreadcnt"]["desc"]) + "\t" +str(statsData["sinkmappingcnt"]["desc"]) + "\t" + str(statsData["totalmappings"]["desc"]) + "\n")
                f_summary.write(source_filename + "\t" + str(statsData["ventreadcnt"]["value"]) + "\t" + str(statsData["mappedreadcnt"]["value"]) + "\t" +str(statsData["sinkmappingcnt"]["value"]) + "\t" + str(statsData["totalmappings"]["value"]) + "\n")
                f_summary.close()
                updateMessages(msgHistory, msgScreen, f_log, logLevel, "FINISHED. Terminating all mapping jobs and ventilator.")
                
                if args.single_host:
                    # stop the mappers
                    print job_ids
                    for jid in job_ids["mappers"]:
                        os.kill(int(jid), signal.SIGINT)
                    # stop the vent
                    os.kill(int(job_ids["vent"]), signal.SIGINT)
                    # wait for sink to finish
                    os.waitpid(int(job_ids["sink"]), 0)
                else:
                    o,e = deleteQSub(job_ids["mappers"])
                    o,e = deleteQSub(job_ids["vent"])               
                if gui:
                    curses.endwin()
                
                time.sleep(10)
                
                sys.exit(0)
        
        elif sender == 'ERROR':
            quitController("ERROR [" + sampleID + "] " +  message, logDirectory, job_ids, f_log, f_summary,f_univerror, f_contigs, gui, args.single_host)
        elif sender == 'WARNING':
            updateMessages(msgHistory, msgScreen, f_log, logLevel, "[WARNING]" + message)
            socket.send("ok")
        else:
            print "unkown sender", message
