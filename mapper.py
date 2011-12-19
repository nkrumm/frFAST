import numpy as np
import operator

class Mappers:
	def __init__(self, numMappers):
		self.numMappers = numMappers
		self.mappers = {}
		self.curr_mapperID = 0
	
	
	def registerMapper(self, node='NA'):
		self.curr_mapperID +=1
		mapperID = self.curr_mapperID	
		self.mappers[int(mapperID)] = {"mapperID": mapperID, \
							 "node": 	 node, \
							 "working":  False, \
							 "rank": 	 0,\
							 "task": ""}
		return mapperID
	
	def startMapper(self, mapperID,seqListSize, mappedSeqCnt):
		self.mappers[int(mapperID)]["seqListSize"] = seqListSize
		self.mappers[int(mapperID)]["mappedSeqCnt"] = mappedSeqCnt
		self.mappers[int(mapperID)]["working"] = True
		self.mappers[int(mapperID)]["rank"] += 1
	
	def updateMapperTask(self, mapperID, task):
		self.mappers[int(mapperID)]["task"] = task

	def updateMapperWorkingStatus(self, mapperID, working):
		self.mappers[int(mapperID)]["working"] = working
	
	
	def updateMapper(self, mapperID, contig, mappingCnt, mappedSeqCnt, task=None, working = None):
		self.mappers[int(mapperID)]["contig"] = contig
		self.mappers[int(mapperID)]["mappingCnt"] = mappingCnt
		self.mappers[int(mapperID)]["mappedSeqCnt"] = mappedSeqCnt
		if task != None:
			self.mappers[int(mapperID)]["task"] = task
		
	def stopMapper(self, mapperID, mappingCnt=0, mappedSeqCnt=0):
		#chunks.append([mapperID, mappers[int(mapperID)]["node"], mappingCnt, mappedSeqCnt])
		del self.mappers[int(mapperID)]
	
	def printMappers(self):
		print "mapperID\tnode\tstatus\ttask\trank"
		for _mapper in self.mappers.values():
			wstr = "Working" if _mapper["working"] else "Free"
			print str(_mapper["mapperID"]) + "\t\t" + str(_mapper["node"]) + "\t" + wstr + "\t" + str(_mapper["task"]) + "\t" + str(_mapper["rank"])
	
	def getNumWorking(self):
		working = np.array(map(operator.itemgetter("working"),self.mappers.values()))
		return np.sum(working)
	
	def getNumRegistered(self):
		return len(self.mappers)
	
	def getNextDestination(self):
		if len(self.mappers) == 0:
			return 0, None # There are no mappers registered yet!
		else:
			#ranks = np.array(map(operator.itemgetter("rank"),self.mappers.values()))
			working = np.array(map(operator.itemgetter("working"),self.mappers.values()))
			mapperIDs = np.array(map(operator.itemgetter("mapperID"),self.mappers.values()))
			if np.sum(working) == len(self.mappers):
				return 0, None # they are all busy
			#print mappers
			available = working==False
			nextMapperID = np.min(mapperIDs[available])
			#print "returning next available mapper to vent: " ,self.mappers[nextMapperID]["node"] + " with mapperID = " + str(nextMapperID)
			return self.mappers[nextMapperID]["node"], nextMapperID

# 
# m = Mappers(3)
# m.registerMapper("hello")
# m.registerMapper("hello2")
# m.registerMapper("hello3")
# m.printMappers()
# m.startMapper(1,100,200)
# m.startMapper(2,100,200)
# m.startMapper(3,100,200)
# m.printMappers()
# m.registerMapper("hello4")
# 
# m.stopMapper(3,100,200)
# m.printMappers()
# 
# m.getNextDestination()
# m.startMapper(4,100,200)
