CXX = g++
CXXFLAGS ?= -ansi  -Winline -Wshadow -Wall -O3
SOURCES = mgsolve.cpp tim.cpp lina.cpp michael.cpp
OBJECTS = mgsolve.o tim.o 
HEADERS = header.h
LIBS = -lm

.PHONY: all clean

all:  mgsolve.o mgsolve

# generic compilation rule
%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $<

#how to link
mgsolve: ${OBJECTS}
	${CXX} -o $@ ${OBJECTS} ${LIBS}

clean:
	rm -f *.o *~ 
