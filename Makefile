CXX      ?= g++
CXXFLAGS ?= -std=c++20
SRCS = src/main.cpp 
HEADERS = src/Matrix.hpp
OBJS = main.o
EXEC = main

.PHONY = all clean

all: ${EXEC}

${EXEC}: ${OBJS}
	${CXX} ${CXXFLAGS} ${OBJS} -o ${EXEC}

${OBJS} : ${SRCS} ${HEADERS}
	${CXX} ${CXXFLAGS} -c ${SRCS}

clean:
	${RM} *.o

distclean: clean
	$(RM) $(EXEC)