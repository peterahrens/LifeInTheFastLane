#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <mpi.h>

#define TIMEOUT 0.1

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

double wall_time (void) {
  struct timeval t;
  gettimeofday(&t, NULL);
  return 1.0*t.tv_sec + 1.0e-6*t.tv_usec;
}

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
  int rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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

  unsigned *initial = NULL;
  if (rank == 0) {
    initial = (unsigned*)malloc(sizeof(unsigned) * height * width);

    for (unsigned y = 0; y < height; y++) {
      for (unsigned x = 0; x < width; x++) {
        initial[y * width + x] = (rand() % 3 == 0);
      }
    }
  }

  double reference_time = 0;
  unsigned *reference;

  if (rank == 0) {
    if (check) {
      for (trials = 1; reference_time < TIMEOUT; trials *= 2) {
        reference = reference_life(height, width, initial, 1);
        free(reference);

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
  }

  double test_time = 0;
  unsigned *test;

  test_time = -wall_time();
  test = life(height, width, initial, iters);
  test_time += wall_time();

  if (rank == 0) {
    if (check) {
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
  }

  MPI_Finalize();

  return 0;
}
