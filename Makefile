CC=         g++
CFLAGS=     -O2

all: main.cpp
	$(CC) $(CFLAGS) main.cpp -o kmers

clean:
	rm kmers
