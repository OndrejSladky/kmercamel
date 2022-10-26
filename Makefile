CC=         g++
CFLAGS=     -std=c++17 -O2

all: main.cpp
	$(CC) $(CFLAGS) main.cpp -o kmers

clean:
	rm -f kmers
