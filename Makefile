.PHONY: all clean test cpptest converttest verify quick-verify

CXX=         g++
CXXFLAGS=    -g -Wall -Wno-unused-function -std=c++17 -O2
LDFLAGS=     -lz
SRC=         src
TESTS=			 tests
GTEST=       $(TESTS)/googletest/googletest

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

kmercamel: $(SRC)/main.cpp $(SRC)/$(wildcard *.cpp *.h *.hpp) src/version.h
	./create-version.sh
	$(CXX) $(CXXFLAGS) $(SRC)/main.cpp -o $@ $(LDFLAGS)
	cp kmercamel  🐫 || true

kmercameltest: $(TESTS)/unittest.cpp gtest-all.o $(SRC)/$(wildcard *.cpp *.h *.hpp) $(TESTS)/$(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include $(TESTS)/unittest.cpp gtest-all.o -pthread -o $@ $(LDFLAGS)

gtest-all.o: $(GTEST)/src/gtest-all.cc $(wildcard *.cpp *.h *.hpp)
	$(CXX) $(CXXFLAGS) -isystem $(GTEST)/include -I $(GTEST)/include -I $(GTEST) -DGTEST_CREATE_SHARED_LIBRARY=1 -c -pthread $(GTEST)/src/gtest-all.cc -o $@

src/version.h: src/version
	./create-version.sh


clean:
	rm -f kmercamel
	rm -f 🐫 || true
	rm -f kmercameltest
	rm -r -f ./bin
	rm -f gtest-all.o
	rm -f src/version.h
