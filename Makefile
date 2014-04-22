CXX ?= g++
CXXFLAGS ?= -ansi  -Winline -Wshadow -O3
SOURCES = tim.cpp lina.cpp michael.cpp
OBJECTS = mgsolve.o tim.o lina.o michael.o
HEADERS = header.h


.PHONY: all clean

all: tim.o lina.o michael.o mgsolve.o mgsolve

#generic compilation rule
%.o : %.cpp
	${CXX} ${CXXFLAGS} -c $<

#how to link
mgsolve: ${OBJECTS}
	${CXX} -o $@ ${OBJECTS} ${LIBS}
	rm -f ${OBJECTS}
	
clean:
	rm -f *.o *~

