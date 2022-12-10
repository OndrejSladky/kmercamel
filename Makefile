.PHONY: all clean test cpptest converttest verify build-win

CXX=         g++
CXXFLAGS=    -g -Wall -Wno-unused-function -std=c++17 -O2
GTEST=       googletest/googletest

all: kmers

# Build kmers so that it can be executed on windows.
build-win: kmers.exe

test: cpptest converttest verify

verify: verify.py kmers
	./verify.py ./spneumoniae.fa

quick-verify: verify.py kmers
	./verify.py --quick ./spneumoniae.fa

cpptest: kmerstest
	./kmerstest

# Execute unittest on Windows.
cpptest-win: kmerstest.exe
	start ./kmerstest.exe

converttest: convert_superstring_unittest.py
	./convert_superstring_unittest.py

kmers: main.cpp $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) main.cpp -o $@

kmers.exe: main.cpp $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) main.cpp -o $@

kmerstest: unittest.cpp gtest-all.o $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include unittest.cpp gtest-all.o -pthread -o $@

kmerstest.exe: unittest.cpp gtest-all.o $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include unittest.cpp gtest-all.o -pthread -o ./$@

gtest-all.o: $(GTEST)/src/gtest-all.cc $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include -I $(GTEST) -DGTEST_CREATE_SHARED_LIBRARY=1 -c -pthread $(GTEST)/src/gtest-all.cc -o $@

clean:
	rm -f kmers
	rm -f kmerstest
	rm -f kmers.exe
	rm -f kmerstest.exe
	rm -r -f ./bin
	rm -f gtest-all.o
