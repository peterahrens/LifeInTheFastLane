CFLAGS = -std=c99 -O3 -march=native

all: simple padded vectorized unrolled
clean:
	rm simple padded vectorized unrolled *.o

reference.o: reference.c
	gcc ${CFLAGS} -c reference.c

simple.o: simple.c
	gcc ${CFLAGS} -c simple.c

padded.o: padded.c
	gcc ${CFLAGS} -c padded.c

vectorized.o: vectorized.c
	gcc ${CFLAGS} -c vectorized.c

unrolled.o: unrolled.c
	gcc ${CFLAGS} -c unrolled.c

simple: bench.c simple.o reference.o
	gcc ${CFLAGS} bench.c reference.o simple.o -o simple

padded: bench.c padded.o reference.o
	gcc ${CFLAGS} bench.c reference.o padded.o -o padded

vectorized: bench.c vectorized.o reference.o
	gcc ${CFLAGS} bench.c reference.o vectorized.o -o vectorized

unrolled: bench.c vectorized.o reference.o unrolled.o
	gcc ${CFLAGS} bench.c reference.o unrolled.o -o unrolled
