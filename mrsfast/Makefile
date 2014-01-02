CC=gcc
CFLAGS = -c -O3 -g
#LDFLAGS = -static -lz -lm
LDPATH = -L/net/gs/vol3/software/modules-sw/zeromq/2.1.11/Linux/RHEL6/x86_64/lib/ 
LDFLAGS = -lz -lm -lzmq
SOURCES = baseFAST.c CommandLineParser.c Common.c HashTable.c MrsFAST.c Output.c Reads.c RefGenome.c 
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = mrsfast

GCC_VERSION = $(shell $(CC) -v 2>&1 | grep version | cut -d' ' -f3  | cut -d'.' -f1)

ifeq ($(GCC_VERSION), 3)
	CC = gcc4
endif


all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o
		
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDPATH) $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@
