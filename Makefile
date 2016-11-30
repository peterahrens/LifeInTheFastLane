# Copyright (c) 2016, Los Alamos National Security, LLC
# All rights reserved.
# Copyright 2016. Los Alamos National Security, LLC. This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.
# Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 1.      Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 2.      Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 3.      Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

CC = gcc
MPICC = mpicc
CFLAGS = -std=c99 -march=native -O3 -lm
OMPFLAGS = -fopenmp

all: simple padded blocked sse2 avx2 streaming omp mpi
clean:
	rm -rf simple padded blocked sse2 avx2 streaming omp mpi *.o *.dSYM

reference.o: reference.c
	${CC} ${CFLAGS} -c reference.c

simple.o: simple.c
	${CC} ${CFLAGS} -c simple.c

padded.o: padded.c
	${CC} ${CFLAGS} -c padded.c

blocked.o: blocked.c
	${CC} ${CFLAGS} -c blocked.c

sse2.o: sse2.c
	${CC} ${CFLAGS} -c sse2.c

avx2.o: avx2.c
	${CC} ${CFLAGS} -c avx2.c

streaming.o: streaming.c
	${CC} ${CFLAGS} -c streaming.c

omp.o: omp.c
	${CC} ${OMPFLAGS} ${CFLAGS} -c omp.c

mpi.o: mpi.c
	${MPICC} ${OMPFLAGS} ${CFLAGS} -c mpi.c

simple: bench.c simple.o reference.o
	${CC} ${CFLAGS} bench.c reference.o simple.o -o simple

padded: bench.c padded.o reference.o
	${CC} ${CFLAGS} bench.c reference.o padded.o -o padded

blocked: bench.c blocked.o reference.o
	${CC} ${CFLAGS} bench.c reference.o blocked.o -o blocked

sse2: bench.c sse2.o reference.o
	${CC} ${CFLAGS} bench.c reference.o sse2.o -o sse2

avx2: bench.c avx2.o reference.o
	${CC} ${CFLAGS} bench.c reference.o avx2.o -o avx2

streaming: bench.c streaming.o reference.o
	${CC} ${CFLAGS} bench.c reference.o streaming.o -o streaming

omp: bench.c reference.o omp.o
	${CC} ${CFLAGS} ${OMPFLAGS} bench.c reference.o omp.o -o omp

mpi: bench_mpi.c reference.o mpi.o
	${MPICC} ${CFLAGS} ${OMPFLAGS} bench_mpi.c reference.o mpi.o -o mpi
