  CXX=g++
  LD=g++ 
  USER_BLAS=gslcblas
  CXXFLAGS=-Wall -pedantic -std=c++11 -O3 -D_GNU_SOURCE_ -g  -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DDBGLEVEL=0 -DPROGNFO -I -Lsrc -DPCGRNG
  LDFLAGS=-L/usr/local/lib -lm -lpthread -lgsl
  # if armadillo
  #LDFLAGS += -larmadillo
  LDFLAGS += -l${USER_BLAS}
  
  CTAGS=ctags > tags
  LIBS=

  METSEGOBJ = src/manopt.o\
	src/memory.o\
	src/debug.o\
	src/info.o\
	src/stringutils.o\
	src/fileio.o\
	src/mm2.o\
	src/mathematics.o\
	src/sort.o\
	src/vstack.o\
	src/segmentstack.o\
	src/mtc.o\
	src/metseg.o
	

all: pmat 

pmat: ${METSEGOBJ}
	${CXX} $(CXXFLAGS) ${METSEGOBJ} -o pmat $(LDFLAGS)

clean:
	rm -rf src/*.o pmat
	rm -f test.mr.DMR

test1: 
	./pmat -t 32 -a A -b B -X 104 -Y 104 -m 5 -d 0.01 ../examples/test.dat > test.mr.DMR


test: test1

.PHONY: all clean
