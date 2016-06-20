all: simple padded vectorized unrolled
clean:
	rm simple padded vectorized unrolled *.o

reference.o: reference.c
	gcc -c reference.c

simple.o: simple.c
	gcc -O3 -c simple.c

padded.o: padded.c
	gcc -O3 -c padded.c

vectorized.o: vectorized.c
	gcc -march=native -O3 -c vectorized.c

unrolled.o: unrolled.c
	gcc -march=native -O3 -c unrolled.c

simple: bench.c simple.o reference.o
	gcc bench.c reference.o simple.o -o simple

padded: bench.c padded.o reference.o
	gcc bench.c reference.o padded.o -o padded

vectorized: bench.c vectorized.o reference.o
	gcc bench.c reference.o vectorized.o -o vectorized

unrolled: bench.c vectorized.o reference.o unrolled.o
	gcc bench.c reference.o unrolled.o -o unrolled
