.PHONY: all clean test cpptest converttest verify

CXX=         g++
CXXFLAGS=    -g -Wall -Wno-unused-function -std=c++17 -O2
GTEST=       googletest/googletest

all: kmers

test: cpptest converttest verify

verify: verify.py kmers
	./verify.py ./spneumoniae.fa

cpptest: kmerstest
	./kmerstest

converttest: convert_superstring_unittest.py
	./convert_superstring_unittest.py

kmers: main.cpp $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) main.cpp -o $@

kmerstest: unittest.cpp gtest-all.o $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -I $(GTEST)/include unittest.cpp gtest-all.o -o $@

gtest-all.o: $(GTEST)/src/gtest-all.cc $(wildcard *.cpp *.h *.hpp)
	git submodule init
	git submodule update
	$(CXX) $(CXXFLAGS)  -I $(GTEST)/include -I $(GTEST) -DGTEST_CREATE_SHARED_LIBRARY=1 -c $(GTEST)/src/gtest-all.cc -o $@

clean:
	rm -f kmers
	rm -f kmerstest
	rm -f gtest-all.o
