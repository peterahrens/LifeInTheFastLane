CFLAGS = -std=c99 -O3 -march=native

all: simple padded sse2 avx2 unrolled
clean:
	rm simple padded sse2 avx2 unrolled *.o

reference.o: reference.c
	gcc ${CFLAGS} -c reference.c

simple.o: simple.c
	gcc ${CFLAGS} -c simple.c

padded.o: padded.c
	gcc ${CFLAGS} -c padded.c

sse2.o: sse2.c
	gcc ${CFLAGS} -c sse2.c

avx2.o: avx2.c
	gcc ${CFLAGS} -c avx2.c

unrolled.o: unrolled.c
	gcc ${CFLAGS} -c unrolled.c

simple: bench.c simple.o reference.o
	gcc ${CFLAGS} bench.c reference.o simple.o -o simple

padded: bench.c padded.o reference.o
	gcc ${CFLAGS} bench.c reference.o padded.o -o padded

sse2: bench.c sse2.o reference.o
	gcc ${CFLAGS} bench.c reference.o sse2.o -o sse2

avx2: bench.c avx2.o reference.o
	gcc ${CFLAGS} bench.c reference.o avx2.o -o avx2

unrolled: bench.c reference.o unrolled.o
	gcc ${CFLAGS} bench.c reference.o unrolled.o -o unrolled
