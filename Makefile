EXTRA_CXXFLAGS=
CXXFLAGS=-O3 -Wall -std=c++17 $(EXTRA_CXXFLAGS)

all: pmpeg.o compress.o decompress.o
	g++ -o pmpeg pmpeg.o compress.o decompress.o

pmpeg.o: pmpeg.cpp
	g++ $(CXXFLAGS) -c pmpeg.cpp

compress.o: compress.cpp
	g++ $(CXXFLAGS) -c compress.cpp

decompress.o: decompress.cpp
	g++ $(CXXFLAGS) -c decompress.cpp

clean:
	rm -f pmpeg *.o
