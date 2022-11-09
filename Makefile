.PHONY: all clean

CXX=         g++
CXXFLAGS=    -std=c++17 -O2

all: kmers

kmers: main.cpp $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) main.cpp -o $@

clean:
	rm -f kmers
