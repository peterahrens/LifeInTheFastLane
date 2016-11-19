/*
 * Copyright (c) 2016, Peter Ahrens
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Los Alamos National Laboratories nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL PETER AHRENS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * LA-CC-17-018
 */

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#define TIMEOUT 0.1

//Show the difference between test results and reference results by highlighting
//places where the two disagree.
static void show_diff (unsigned height, unsigned width, unsigned *universe, unsigned *ref) {
  printf("+");
  for (int x = 0; x < width; x++) {
    printf("--");
  }
  printf("+\n");
  for (int y = 0; y < height; y++) {
    printf("|");
    for (int x = 0; x < width; x++) {
      printf(universe[y * width + x] != ref[y * width + x] ? "##" : "  ");
    }
    printf("|\n");
  }
  printf("+");
  for (int x = 0; x < width; x++) {
    printf("--");
  }
  printf("+\n");
}

//Show a particular life grid.
static void show (unsigned height, unsigned width, unsigned *universe) {
  printf("+");
  for (int x = 0; x < width; x++) {
      printf("--");
  }
  printf("+\n");
  for (int y = 0; y < height; y++) {
    printf("|");
    for (int x = 0; x < width; x++) {
      printf(universe[y * width + x] ? "##" : "  ");
    }
    printf("|\n");
  }
  printf("+");
  for (int x = 0; x < width; x++) {
      printf("--");
  }
  printf("+\n");
}

//Return time time of day as a double-precision floating point value.
double wall_time (void) {
  struct timeval t;
  //It is critically important to use an accurate timer. Many common functions
  //that return the time are not accurate enough for timing code. Since timers
  //are typically system-specific, research timers for your system.
  //Surprisingly, omp_get_wtime() is quite good and is available everywhere
  //there is OpenMP.
  gettimeofday(&t, NULL);
  return 1.0*t.tv_sec + 1.0e-6*t.tv_usec;
}

//We will link in different methods for life.
unsigned *life (const unsigned height,
                const unsigned width,
                const unsigned * const initial,
                const unsigned iters);

unsigned *reference_life (unsigned height,
                          unsigned width,
                          unsigned *initial,
                          unsigned iters);

//Calling conventions for bench.c executables:
//The first number (if included) is the height
//The second number (if included) is the width
//The third number (if included) is the number of iterations
//The fourth number (if included) is 1 if you want to compare to the
//  reference implementation and 0 (default) if you do not. This option
//  is included because the reference implementation can be very slow.
int main (int argc, char **argv) {
  unsigned height = 0;
  unsigned width = 0;
  unsigned iters = 0;
  unsigned check = 0;
  unsigned trials;

  if (argc > 2) {
    width = atoi(argv[1]);
    height = atoi(argv[2]);
  }
  if (argc > 3) {
    iters = atoi(argv[3]);
  }
  if (argc > 4) {
    check = atoi(argv[3]);
  }

  //Default values.
  if (width <= 0) {
    width = 8192;
  }
  if (height <= 0) {
    height = 8192;
  }
  if (iters <= 0) {
    iters = 256;
  }
  if (check <= 0) {
    check = 0;
  }

  unsigned *initial = (unsigned*)malloc(sizeof(unsigned) * height * width);

  //Initialize universe randomly.
  for (unsigned y = 0; y < height; y++) {
    for (unsigned x = 0; x < width; x++) {
      initial[y * width + x] = (rand() % 3 == 0);
    }
  }

  double reference_time = 0;
  unsigned *reference;
  if (check) {

    //We must run the benchmarking application for a sufficient length of time
    //to avoid small variations in processing speed. We do this by running an
    //increasing number of trials until it takes at least TIMEOUT seconds.
    for (trials = 1; reference_time < TIMEOUT; trials *= 2) {

      //Unless you want to measure the cache warm-up time, it is usually a good
      //idea to run the problem for one iteration first to load the problem
      //into cache.
      reference = reference_life(height, width, initial, 1);
      free(reference);

      //Benchmark "trials" runs of life.
      reference_time = -wall_time();
      for (int i = 0; i < trials - 1; ++i){
        reference = reference_life(height, width, initial, iters);
        free(reference);
      }
      reference = reference_life(height, width, initial, iters);
      reference_time += wall_time();
    }
    trials /= 2;
    reference_time /= trials;
  }

  double test_time = 0;
  unsigned *test;

  //We must run the benchmarking application for a sufficient length of time
  //to avoid small variations in processing speed. We do this by running an
  //increasing number of trials until it takes at least TIMEOUT seconds.
  for (trials = 1; test_time < TIMEOUT; trials *= 2) {

    //Unless you want to measure the cache warm-up time, it is usually a good
    //idea to run the problem for one iteration first to load the problem
    //into cache.
    test = life(height, width, initial, 1);
    free(test);

    //Benchmark "trials" runs of life.
    test_time = -wall_time();
    for (int i = 0; i < trials - 1; ++i){
      test = life(height, width, initial, iters);
      free(test);
    }
    test = life(height, width, initial, iters);
    test_time += wall_time();
  }
  trials /= 2;
  test_time /= trials;

  if (check) {
    //Ensure that reference and test results are identical.
    for (unsigned y = 0; y < height; y++) {
      for (unsigned x = 0; x < width; x++) {
        if (test[y * width + x] != reference[y * width + x]) {
          printf("Test life results do not match reference.\n");
          show(height, width, test);
          show(height, width, reference);
          show_diff(height, width, test, reference);
          exit(-1);
        }
      }
    }
  }

  if (check) {
    printf("Reference time: %g\n", reference_time);
  }
  printf("     Test time: %g\n", test_time);
  if (check) {
    printf("       Speedup: %g\n", reference_time / test_time);
  }

  free(initial);
  if (check) {
    free(reference);
  }
  free(test);

  return 0;
}
