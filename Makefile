CC = g++
FLAGS = -std=c++11 -fPIC -Wall
LFLAGS = -L. -lGFW

all: libGFW.so Test
Test: libGFW.so
	$(CC) $(FLAGS) $(LFLAGS) -o Test Test.C
libGFW.so: GFWCumulant.o GFW.o
	$(CC) $(FLAGS) -shared -o libGFW.so GFW.o GFWCumulant.o
GFWCumulant.o:
	$(CC) $(FLAGS) -c -o GFWCumulant.o GFWCumulant.cxx
GFW.o:
	$(CC) $(FLAGS) -c -o GFW.o GFW.cxx
clean:
	rm *.o *.so Test
