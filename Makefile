CC = gcc
CFLAGS = -std=c99 -march=native -O3
OMPFLAGS = -fopenmp

all: simple padded sse2 avx2 omp
clean:
	rm simple padded sse2 avx2 omp *.o

reference.o: reference.c
	${CC} ${CFLAGS} -c reference.c

simple.o: simple.c
	${CC} ${CFLAGS} -c simple.c

padded.o: padded.c
	${CC} ${CFLAGS} -c padded.c

sse2.o: sse2.c
	${CC} ${CFLAGS} -c sse2.c

avx2.o: avx2.c
	${CC} ${CFLAGS} -c avx2.c

omp.o: omp.c
	${CC} ${OMPFLAGS} ${CFLAGS} -c omp.c

simple: bench.c simple.o reference.o
	${CC} ${CFLAGS} bench.c reference.o simple.o -o simple

padded: bench.c padded.o reference.o
	${CC} ${CFLAGS} bench.c reference.o padded.o -o padded

sse2: bench.c sse2.o reference.o
	${CC} ${CFLAGS} bench.c reference.o sse2.o -o sse2

avx2: bench.c avx2.o reference.o
	${CC} ${CFLAGS} bench.c reference.o avx2.o -o avx2

omp: bench.c reference.o omp.o
	${CC} ${CFLAGS} ${OMPFLAGS} bench.c reference.o omp.o -o omp
