CC = g++
FLAGS = -std=c++11 -fPIC -Wall
LFLAGS = -L. -lGFW

all: libGFW.so Test
Test: libGFW.so Test.C
	$(CC) $(FLAGS) $(LFLAGS) -o Test Test.C
libGFW.so: GFWCumulant.o GFW.o
	$(CC) $(FLAGS) -shared -o libGFW.so GFW.o GFWCumulant.o
GFWCumulant.o: GFWCumulant.cxx GFWCumulant.h
	$(CC) $(FLAGS) -c -o GFWCumulant.o GFWCumulant.cxx
GFW.o: GFW.cxx GFW.h
	$(CC) $(FLAGS) -c -o GFW.o GFW.cxx
clean:
	rm *.o *.so Test
