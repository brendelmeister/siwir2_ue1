CXX = g++
CXXFLAGS ?= -ansi -Winline -Wshadow -Wall -O3 -ggdb
SOURCES = mgsolveN.cpp
OBJECTS = mgsolveN.o
HEADERS = header.h
LIBS = -lm

.PHONY: all clean

all: mgsolveN.o mgsolveN

# generic compilation rule
%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $<

#how to link
mgsolveN: ${OBJECTS}
	${CXX} -o $@ ${OBJECTS} ${LIBS}

clean:
	rm -f *.o *~ *.txt mgsolveN

