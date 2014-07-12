CC 		= g++
CCFLAGS = -Wall -O3 -fopenmp #-lrt
IFLAGS 	= -I/opt/local/include 
LFLAGS 	= -lgmpxx -lgmp

all: Multimap.cpp Multimap.h main.cpp Makefile
	$(CC) $(CCFLAGS) $(IFLAGS) -c -o Multimap.o Multimap.cpp
	$(CC) $(CCFLAGS) $(IFLAGS) $(LFLAGS) main.cpp Multimap.o -o multimap

clean:
	rm -f *.o
	rm -f multimap
