.PHONY: all clean

CC=         g++
CFLAGS=     -std=c++17 -O2

all: kmers

kmers: main.cpp $(wildcard *.cpp *.h *.hpp)
	$(CC) $(CFLAGS) main.cpp -o $@

clean:
	rm -f kmers
