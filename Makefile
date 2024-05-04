CXX      ?= g++
CXXFLAGS ?= -std=c++20 -O3
SRCS = src/main.cpp 
HEADERS = src/Matrix.hpp
OBJS = main.o
EXEC = main
PACS_ROOT = ../../pacs-examples/Examples
UTILITIES = ${PACS_ROOT}/src/Utilities

.PHONY = all clean

all: ${EXEC}

${EXEC}: ${OBJS}
	${CXX} ${CXXFLAGS} ${OBJS} -o ${EXEC}

${OBJS} : ${SRCS} ${HEADERS}
	${CXX} ${CXXFLAGS} -c -I${UTILITIES} ${SRCS}

clean:
	${RM} *.o

distclean: clean
	$(RM) $(EXEC)