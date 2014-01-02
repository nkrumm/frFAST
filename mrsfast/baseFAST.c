/*
 * Copyright (c) <2008 - 2009>, University of Washington, Simon Fraser University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the name of the <ORGANIZATION> nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Author         : Faraz Hach
 * Email          : fhach AT cs DOT sfu
 * Last Update    : 2009-02-01
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"
#include "CommandLineParser.h"
#include "Reads.h"
#include "Output.h"
#include "HashTable.h"
#include "MrsFAST.h"
#include <zmq.h>
#include "zhelpers.h"
#include <time.h>
#include <errno.h>
#include <assert.h>

void *context;
void *requester; 
char *mapperID;
char *sinkNode;
char *ventNode;
void *out_sock;
void *subsocket;

char 				*versionNumber = "2.3";			// Current Version
unsigned char		seqFastq;



int doMrsFAST(){
		//configHashTable();
		

		Read *seqList;
		unsigned int seqListSize;
		int fc;
		int samplingLocsSize;
		int *samplingLocs;
		double totalLoadingTime = 0;
		double totalMappingTime = 0;
		double startTime;
		double loadingTime;
		double mappingTime;
		double lstartTime;
		double ppTime;
		double tmpTime;
		char *prevGen = getMem(CONTIG_NAME_SIZE);
		prevGen[0]='\0';
		char *curGen;
		int	flag;
		double maxMem=0;
		char fname1[FILE_NAME_LENGTH];
		char fname2[FILE_NAME_LENGTH];
		char fname3[FILE_NAME_LENGTH];
		char fname4[FILE_NAME_LENGTH];
		char fname5[FILE_NAME_LENGTH];
		
		
		clearHashes();
		
		// used for ZMQ communication
		char out_msg [255];
		int n;	
		
		
		// Loading Sequences & Sampling Locations
		startTime = getTime();
		
		if (!readAllReads(seqFile1, seqFile2, seqCompressed, &seqFastq, pairedEndMode, &seqList, &seqListSize))
			{
				return 1;
			}
		
		//loadSamplingLocations(&samplingLocs, &samplingLocsSize);
		totalLoadingTime += getTime()-startTime;


		// Preparing output
		initOutput(mappingOutput, outCompressed);

		fprintf(stdout, "-----------------------------------------------------------------------------------------------------------\n");
		fprintf(stdout, "| %15s | %15s | %15s | %15s | %15s %15s |\n","Genome Name","Loading Time", "Mapping Time", "Memory Usage(M)","Total Mappings","Mapped reads");
		fprintf(stdout, "-----------------------------------------------------------------------------------------------------------\n");
		
		//fprintf(stdout, "MAPPER START %s %d %d\n",mapperID, seqListSize, mappedSeqCnt);

		sprintf (out_msg, "MAPPER START %s %d %d",mapperID, seqListSize, mappedSeqCnt);
		s_send (requester, out_msg);
		s_recv (requester);
		
		
		for (fc = 0; fc < fileCnt; fc++)
		{
			if (!initLoadingHashTable(fileName[fc][1]))
			{
				return 1;
			}
			loadSamplingLocations(&samplingLocs, &samplingLocsSize);
			mappingTime = 0;
			loadingTime = 0;
			prevGen[0] = '\0';
			flag = 1;
			do
			{

				flag = loadHashTable ( &tmpTime );  			// Reading a fragment
				
				curGen = getRefGenomeName();

				
				// First Time
				if (flag && prevGen[0]== '\0')
				{
					sprintf(prevGen, "%s", curGen);

				}

				if ( !flag || strcmp(prevGen, curGen)!=0)
				{
					fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n",
							prevGen,loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
					fflush(stdout);
					
					sprintf (out_msg, "MAPPER UPDATE %s %s %d %d",mapperID, prevGen, mappingCnt, mappedSeqCnt);
					s_send (requester, out_msg);
					s_recv (requester);
					
					totalMappingTime += mappingTime;
					totalLoadingTime += loadingTime;

					loadingTime = 0;
					mappingTime = 0;
					maxMem = 0;

					if (!flag)
					{
						break;
					}
				}
				else if (progressRep && mappingTime != 0)
				{
					fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n",
							prevGen,loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
					fflush(stdout);
				}


				sprintf(prevGen, "%s", curGen);

				loadingTime += tmpTime;
				lstartTime = getTime();


				initFAST(seqList, seqListSize, samplingLocs, samplingLocsSize, fileName[fc][0]);
				
				mapSingleEndSeq();

				mappingTime += getTime() - lstartTime;
				if (maxMem < getMemUsage())
				{
					maxMem = getMemUsage();
				}

			} while (flag);

		} // end for;
		
		finalizeFAST();
		finalizeLoadingHashTable();


		finalizeOutput();

		fprintf(stdout, "-----------------------------------------------------------------------------------------------------------\n");
		fprintf(stdout, "%19s%16.2f%18.2f\n\n", "Total:",totalLoadingTime, totalMappingTime);
		if (pairedEndDiscordantMode)
			fprintf(stdout, "Post Processing Time: %18.2f \n", ppTime);
		fprintf(stdout, "%-30s%10.2f\n","Total Time:", totalMappingTime+totalLoadingTime);
		fprintf(stdout, "%-30s%10d\n","Total No. of Reads:", seqListSize);
		fprintf(stdout, "%-30s%10d\n","Total No. of Mapped Reads:", mappedSeqCnt);
		fprintf(stdout, "%-30s%10lld\n","Total No. of Mappings:", mappingCnt);
		fprintf(stdout, "%-30s%10.0f\n\n","Avg No. of locations verified:", ceil((float)verificationCnt/seqListSize));

		char out_mapping [255];
		
		
		sprintf (out_mapping, "%s %d %d", mapperID,  mappingCnt , mappedSeqCnt);
		s_send (out_sock, out_mapping);
		
		
		
		sprintf (out_msg, "MAPPER FINISH %s %d %d",mapperID, mappingCnt , mappedSeqCnt);
		s_send (requester, out_msg);
		s_recv (requester);
		
		
		
		int cof = (pairedEndMode)?2:1;

		if (progressRep && maxHits != 0)
		{
			int frequency[maxHits+1];
			int i;
			for ( i=0 ; i <= maxHits; i++)
			{
				frequency[i] = 0;
			}


			for (fc = 0; fc < seqListSize; fc++)
			{
				frequency[*(seqList[fc*cof].hits)]++;
			}
			frequency[maxHits] = completedSeqCnt;
			for ( i=0 ; i <= maxHits; i++)
			{
				fprintf(stdout, "%-30s%10d%10d%10.2f%%\n","Reads Mapped to ", i, frequency[i], 100*(float)frequency[i]/(float)seqListSize);
			}
		}

		freeMem(prevGen, CONTIG_NAME_SIZE);
		
		
		//free(seqList);
		//free(samplingLocs);
		//free(prevGen);
		//free(curGen);
		return 1;
}


int main(int argc, char *argv[])
{
	int major, minor, patch;
	zmq_version (&major, &minor, &patch); printf ("ZMQ version is %d.%d.%d\n", major, minor, patch);
	
	if (!parseCommandLine(argc, argv))
		return 1;
	int e;
	context = zmq_init (3);
	errno = 0;
	
	fprintf(stdout, "Controller node: %s\n", controllerNode);
	fprintf(stdout, "Controller port: %d\n", controllerPort);
	
	printf ("Connecting to controller...\n");
 	requester = zmq_socket (context, ZMQ_REQ);
	

	char controlleraddress[255];
	sprintf(controlleraddress, "tcp://%s:%d", controllerNode, controllerPort);
 	fprintf(stdout, "Connecting to controller @ %s\n", controlleraddress); 	
 	e = zmq_connect (requester, controlleraddress);
 	fprintf (stdout, "error %d; error = %s\n", e, zmq_strerror (errno)); 	
 	
 	printf ("Getting mapperID\n");

 	// SEND HOSTNAME 
 	char hostname[100];
	hostname[99] = '\0';
	gethostname(hostname, 1023);
	char message[100] = "MAPPER REGISTER ";
	strcat(message, hostname);
 	e = s_send (requester, message);
 	
 	//fprintf (stdout, "error %d; error = %s\n", e, zmq_strerror (errno)); 	
	
 	// get mapperID
 	mapperID = s_recv (requester);
	fprintf (stdout, "received reply, mapperID = %s \n", mapperID);
	
	e = s_send (requester, "MAPPER GETSINKADDRESS");
 	//fprintf (stdout, "error %d; error = %s\n", e, zmq_strerror (errno));
 	sinkNode = s_recv (requester);
	fprintf (stdout, "got sink address: = %s \n", sinkNode);
 	
 	s_send (requester, "MAPPER GETVENTADDRESS");
 	ventNode = s_recv (requester);
	fprintf (stdout, "got vent address: = %s \n", ventNode);	
	
 	s_send (requester, "MAPPER GETINDEX");
	
	char *indexFileNetwork = NULL;
 	indexFileNetwork = s_recv (requester);
 	fprintf (stdout, "Received index location: %s\n", indexFileNetwork);
	sprintf(fileName[fileCnt][0], "%s",indexFileNetwork);
	sprintf(fileName[fileCnt][1], "%s.index", fileName[fileCnt][0]); 
	fileCnt++;
	
	fprintf(stdout, "Connecting to sink @ %s...\n", sinkNode);
	out_sock = zmq_socket (context, ZMQ_PUSH);
	zmq_connect (out_sock, sinkNode);
	
	
	//fprintf(stdout, "Connecting to controller subscription service @ %s\n", controlleraddress);
	
	char filter[20] = "";
	strcat(filter, mapperID);
	strcat(filter, " ");
	
	
	subsocket = zmq_socket (context, ZMQ_SUB);
	
	char controllersubaddress[255];
	int PUB_SUB_PORT;
	PUB_SUB_PORT = controllerPort + 1; // formerly 7000
	sprintf(controllersubaddress, "tcp://%s:%d", controllerNode, PUB_SUB_PORT);
	
	fprintf(stdout, "Connecting to controller subscription service @ %s\n", controllersubaddress);
 	
	zmq_connect (subsocket, controllersubaddress);
	zmq_setsockopt (subsocket, ZMQ_SUBSCRIBE, filter, strlen(filter));
	
	
	
	char *msg = s_recv (subsocket);
	
	char killsig[255];
	strcpy(killsig,mapperID);
	strcat(killsig," KILL");
	
	if (strcmp(msg,killsig)==1){
		free(msg);
		printf("KILL received from controller. Terminating!\n");
		return 0; // end the while loop and terminate the mapper
	}
	else{
		printf ("%s\n", msg);
		free(msg);
		printf("GO ahead signal for controller!\n");
		int r;
		r = doMrsFAST();
		//close all sockets, zmq will ensure that all unsent messages are sent first
		
		zmq_close(out_sock);
		zmq_close(requester);
		zmq_close(subsocket);
		zmq_term(context);

		return 1;
	}
}
