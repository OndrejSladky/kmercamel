.PHONY: all clean test cpptest converttest verify quick-verify

CXX=         g++
CXXFLAGS=    -g -Wall -Wno-unused-function -std=c++17 -O2
GTEST=       googletest/googletest

all: kmercamel

test: cpptest converttest verify

verify: verify.py kmercamel
	./verify.py ./spneumoniae.fa

quick-verify: verify.py kmercamel
	./verify.py --quick ./spneumoniae.fa

cpptest: kmercameltest
	./kmercameltest

converttest: convert_superstring_unittest.py
	./convert_superstring_unittest.py

kmercamel: main.cpp $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) main.cpp -o $@
	cp kmercamel 🐫 || true

kmercameltest: unittest.cpp gtest-all.o $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include unittest.cpp gtest-all.o -pthread -o $@

gtest-all.o: $(GTEST)/src/gtest-all.cc $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include -I $(GTEST) -DGTEST_CREATE_SHARED_LIBRARY=1 -c -pthread $(GTEST)/src/gtest-all.cc -o $@

clean:
	rm -f kmercamel
	rm -f 🐫 || true
	rm -f kmercameltest
	rm -r -f ./bin
	rm -f gtest-all.o
