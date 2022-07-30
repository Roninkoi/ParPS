CC=mpicc
CCFLAGS=-std=gnu11 -lm

HEADERS = src/util.h src/poisson.h src/solver.h
SOURCES = src/poisson.c src/solver.c src/main.c

.PHONY: clean all
.DEFAULT: all

all: release

clean:
	rm -f parps *.o

release: $(HEADERS) $(SOURCES)
	$(CC) $(CCFLAGS) $(SOURCES) -o parps

