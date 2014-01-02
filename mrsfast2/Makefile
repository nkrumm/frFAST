ALL: mrsfast

#LDFLAGS=-s -static
LDPATH=-L/net/gs/vol3/software/modules-sw/zeromq/2.1.11/Linux/RHEL6/x86_64/lib/
LIBS=-lz -lm -lzmq
CFLAGS= -O2 -s 

mrsfast: baseFAST.o MrsFAST.o Common.o CommandLineParser.o RefGenome.o HashTable.o  Reads.o Output.o
	gcc $^ -o $@ ${LDFLAGS} ${LDPATH} ${LIBS}
	rm *.o
