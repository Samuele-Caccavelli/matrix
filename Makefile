CXX      ?= g++
CXXFLAGS ?= -std=c++20 -O3
SRCS = src/main.cpp 
HEADERS = src/Matrix.hpp
OBJS = main.o
EXEC = main
PACS_ROOT = ../../pacs-examples/Examples
UTILITIES = ${PACS_ROOT}/src/Utilities
MATRIX = lnsp_131.mtx
MATRIX_COMPRESSED = lnsp_131.mtx.gz
MATRIX_LINK = https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/lns/lnsp_131.mtx.gz

.PHONY = all clean clean_matrix

all: ${EXEC}

${EXEC}: ${OBJS} ${MATRIX}
	${CXX} ${CXXFLAGS} ${OBJS} -o ${EXEC}

${OBJS} : ${SRCS} ${HEADERS}
	${CXX} ${CXXFLAGS} -c -I${UTILITIES} ${SRCS}

${MATRIX}: ${MATRIX_COMPRESSED}
	gzip -d ${MATRIX_COMPRESSED}

${MATRIX_COMPRESSED}:
	wget -O ${MATRIX_COMPRESSED} ${MATRIX_LINK}

clean:
	${RM} *.o

clean_matrix:
	${RM} ${MATRIX}

distclean: clean
	$(RM) $(EXEC)