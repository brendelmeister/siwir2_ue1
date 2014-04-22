CXX = g++
CXXFLAGS ?= -ansi  -Winline -Wshadow -O3
SOURCES = mgsolve.cpp tim.cpp lina.cpp michael.cpp
OBJECTS = mgsolve.o tim.o lina.o michael.o
HEADERS = header.h
LIBS = -lm

.PHONY: all clean

all: tim.o lina.o michael.o mgsolve.o mgsolve

# generic compilation rule
%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $<

#how to link
mgsolve: ${OBJECTS}
	${CXX} -o $@ ${OBJECTS} ${LIBS}

clean:
	rm -f *.o *~ 