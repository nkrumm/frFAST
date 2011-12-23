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
 * Last Update    : 2009-12-08
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include "Common.h"
#include "Reads.h"
#include <zmq.h>
#include "zhelpers.h"
#include <time.h>


FILE *_r_fp1;
FILE *_r_fp2;
gzFile _r_gzfp1;
gzFile _r_gzfp2;
Read *_r_seq;
int _r_seqCnt;
int *_r_samplingLocs;

/**********************************************/
char *(*readFirstSeq)(char *);
char *(*readSecondSeq)(char *);
/**********************************************/
char *readFirstSeqTXT( char *seq )
{
	return fgets(seq, SEQ_MAX_LENGTH, _r_fp1);
}

/**********************************************/
char *readSecondSeqTXT( char *seq )
{
	return fgets(seq, SEQ_MAX_LENGTH, _r_fp2);
}
/**********************************************/
char *readFirstSeqGZ( char *seq )
{
	return gzgets(_r_gzfp1, seq, SEQ_MAX_LENGTH);
}

/**********************************************/
char *readSecondSeqGZ( char *seq )
{
	return gzgets(_r_gzfp2, seq, SEQ_MAX_LENGTH);
}
/**********************************************/
int readAllReads(char *fileName1,
						char *fileName2,
						int compressed,
						unsigned char *fastq,
						unsigned char pairedEnd,
						Read **seqList,
						unsigned int *seqListSize)
{
	double startTime=getTime();

	char seq1[SEQ_MAX_LENGTH];
	char rseq1[SEQ_MAX_LENGTH];
	char name1[SEQ_MAX_LENGTH];
	char qual1[SEQ_MAX_LENGTH];
	char seq2[SEQ_MAX_LENGTH];
	char rseq2[SEQ_MAX_LENGTH];
	char name2[SEQ_MAX_LENGTH];
	char qual2[SEQ_MAX_LENGTH];

	char dummy[SEQ_MAX_LENGTH];
	char ch;
	int err1, err2;
	int nCnt;
	int discarded = 0;
	int seqCnt = 0;
	int maxCnt = 0;
	int i;
	Read *list = NULL;

	// connect pair socket
	printf("attempting to connect to ventilator @ %s\n", ventNode);
	
	void *readreceiver = zmq_socket (context, ZMQ_PULL);
    zmq_connect (readreceiver, ventNode);
    char *message = s_recv (readreceiver);
    
    maxCnt = atoi(message);
    printf("... succesfull!\n");    
    printf("ready to receive %d reads\n", maxCnt);
    free (message);

	
	*fastq = 0;
	
	list = getMem(sizeof(Read)*maxCnt);

	// start receiving
	
	int recvdCnt = 0;	
	
	while (1){
		
		char *message = s_recv (readreceiver);
		char *p = NULL;
		
		p = strtok (message," ");
		sprintf(name1, p);
		p = strtok (NULL," ");		
		sprintf(seq1, p);
		//printf(n)
		
		free(message);
		
		if (strcmp(name1,"DONE")==0){ break;} // get out of loop if the last read is received
		
		sprintf(qual1, "*"); // ignore all quality for now
		recvdCnt++;
		
		for (i=0; i<strlen(name1); i++){if (name1[i] == '\n') {name1[i]='\0';}}
		for (i=0; i<strlen(seq1); i++){if (seq1[i] == '\n') {seq1[i]='\0';}}		
		
		//printf("%d: %s %s\n", recvdCnt, name1, seq1);
		
		//count all the Ns in the sequence and set err1 appropriately
		nCnt = 0;
		err1 = 0;
		for (i=0; i<strlen(seq1); i++)
		{
			if (seq1[i] == 'N') {nCnt++;}
		}

		if (nCnt > errThreshold){
			err1 = 1;
		}
		
		
		if (!err1){
			int _mtmp = strlen(seq1);
			list[seqCnt].hits = getMem (1+3*_mtmp+3+strlen(name1)+1);
			list[seqCnt].seq = list[seqCnt].hits + 1;
			list[seqCnt].rseq = list[seqCnt].seq + _mtmp+1;
			list[seqCnt].qual = list[seqCnt].rseq + _mtmp+1;
			list[seqCnt].name = list[seqCnt].qual + _mtmp+1;
	
	
			reverseComplete(seq1, rseq1, _mtmp);
			rseq1[_mtmp] =  '\0';
			int i;
	
			list[seqCnt].hits[0] = 0;
	
			for (i=0; i<=_mtmp; i++)
			{
				list[seqCnt].seq[i] = seq1[i];
				list[seqCnt].rseq[i] = rseq1[i] ;
				list[seqCnt].qual[i] = qual1[i];
			}
			sprintf(list[seqCnt].name,"%s%c", ((char*)name1)+1,'\0');
	
			seqCnt++;
		}
		else { // sequence had too many Ns in it!
			discarded++;
		}
		
// TO DO ADD BREAK CONIDITION FOR LAST PART OF SEQUENCE
// 		if (recvdCnt == maxCnt) { break;} // break out of while(1) loop
	
	}
	
	// CLOSE PAIR SOCKET
	zmq_close (readreceiver);
	
	
	if (seqCnt > 0)
	{
		QUAL_LENGTH = SEQ_LENGTH = strlen(list[0].seq);
		if (! *fastq)
		{
			QUAL_LENGTH = 1;
		}
		//fprintf(stderr, "%d %d\n", SEQ_LENGTH, QUAL_LENGTH);
	}
	else
	{
		fprintf(stdout, "ERR: No reads can be found for mapping\n");
		fprintf(stdout, "recvdCnt: %d\n", recvdCnt);
		fprintf(stdout, "seqCnt: %d\n", seqCnt);
		fprintf(stdout, "discarded: %d\n", discarded);
		char out_msg [255];
		// tell controller that no reads were received!
		sprintf (out_msg, "MAPPER INPUT %s %d %d",mapperID, seqCnt, discarded);
		s_send (requester, out_msg);
		s_recv (requester);
		
		return 0;
	}


	if (pairedEnd)
	{
//		seqCnt /= 2;
	}

/*
//  Closing Files
 	if (!compressed)
 	{
 		fclose(_r_fp1);
 		if ( pairedEnd && fileName2 != NULL )
 		{
 			fclose(_r_fp2);
 		}
 	}
 	else
 	{
 		gzclose(_r_gzfp1);
 		if ( pairedEnd && fileName2 != NULL)
 		{
 			gzclose(_r_fp2);
 		}
 	}
*/
	*seqList = list;
	*seqListSize = seqCnt;

	_r_seq = list;
	_r_seqCnt = seqCnt;
	
	
	
	fprintf(stdout, "%d sequences are read in %0.2f. (%d discarded) [Mem:%0.2f M]\n", seqCnt, (getTime()-startTime), discarded, getMemUsage());
	//totalLoadingTime+=getTime()-startTime;
	
	char out_msg [255];
	sprintf (out_msg, "MAPPER INPUT %s %d %d",mapperID, seqCnt, discarded);
	s_send (requester, out_msg);
	s_recv (requester);
	
	//zmq_close(readreceiver);
	
	return 1;
}
/**********************************************/
void loadSamplingLocations(int **samplingLocs, int * samplingLocsSize)
{
	int i;
	int samLocsSize = errThreshold + 1;
	int *samLocs = getMem(sizeof(int)*samLocsSize);

	for (i=0; i<samLocsSize; i++)
	{
		samLocs[i] = (SEQ_LENGTH / samLocsSize) *i;
		if ( samLocs[i] + WINDOW_SIZE > SEQ_LENGTH)
			samLocs[i] = SEQ_LENGTH - WINDOW_SIZE;
	}

	// Outputing the sampling locations

/*	int j;
 	for (i=0; i<SEQ_LENGTH; i++)
	{
		fprintf(stdout, "-");
	}
	fprintf(stdout, "\n");

	for ( i=0; i<samLocsSize; i++ )
	{
		for ( j=0; j<samLocs[i]; j++ )
		{
			fprintf(stdout," ");
		}
		for (j=0; j<WINDOW_SIZE; j++)
		{
			fprintf(stdout,"+");
		}
		fprintf(stdout, "\n");
		fflush(stdout);
	}
	for ( i=0; i<SEQ_LENGTH; i++ )
	{
		fprintf(stdout, "-");
	}
	fprintf(stdout, "\n");*/
	*samplingLocs = samLocs;
	*samplingLocsSize = samLocsSize;
	_r_samplingLocs = samLocs;
}

void finalizeReads(char *fileName)
{
	FILE *fp1=NULL;

	if (fileName != NULL)
	{
		fp1 = fileOpen(fileName, "w");
	}
	if (pairedEndMode)
		_r_seqCnt /=2;

	int i=0;
	for (i = 0; i < _r_seqCnt; i++)
	{
		if (pairedEndMode && _r_seq[2*i].hits[0] == 0 &&  strcmp(_r_seq[2*i].qual,"*")!=0)
		{
			fprintf(fp1,"@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].qual, _r_seq[i*2].name, _r_seq[i*2+1].seq, _r_seq[i*2+1].qual);
		}
		else if (pairedEndMode && _r_seq[2*i].hits[0] == 0)
		{
			fprintf(fp1, ">%s/1\n%s\n>%s/2\n%s\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].name, _r_seq[i*2+1].seq);
		}
		else if (_r_seq[i].hits[0] == 0 && strcmp(_r_seq[i].qual, "*")!=0)
		{
			fprintf(fp1,"@%s\n%s\n+\n%s\n", _r_seq[i].name, _r_seq[i].seq, _r_seq[i].qual);
		}
		else if (_r_seq[i].hits[0] == 0)
		{
			fprintf(fp1,">%s\n%s\n", _r_seq[i].name, _r_seq[i].seq);
		}
	}

	fclose(fp1);
	if (pairedEndMode)
		_r_seqCnt *= 2;

	for (i = 0; i < _r_seqCnt; i++)
	{
		freeMem(_r_seq[i].hits,0);
	}


	freeMem(_r_seq,0);
	freeMem(_r_samplingLocs,0);
}
