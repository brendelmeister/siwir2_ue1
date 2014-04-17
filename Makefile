CXX ?= g++
CXXFLAGS ?= -ansi -Wall -Winline -Wshadow -O3

.PHONY: all clean 

all: mgsolve
clean:
	rm -f  mgsolve.o mgsolve *.txt
	rm -f test.o test


mgsolve: mgsolve.o
	$(CXX) $(CXXFLAGS) -o $@ $^ 

mgsolve.o: mgsolve.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<


test: test.o
	$(CXX) $(CXXFLAGS) -o $@ $^

test.o: test.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<
