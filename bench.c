#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <omp.h>

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
                const unsigned iters,
                const unsigned display);

unsigned *reference_life (unsigned height,
                          unsigned width,
                          unsigned *initial,
                          unsigned iters,
                          unsigned display);

int main (int argc, char **argv) {
  unsigned height = 0;
  unsigned width = 0;
  unsigned iters = 0;
  unsigned display = 0;
  unsigned trials;

  if (argc > 2) {
    width = atoi(argv[1]);
    height = atoi(argv[2]);
  }
  if (argc > 3) {
    iters = atoi(argv[3]);
  }
  if (argc > 4) {
    display = atoi(argv[4]);
  }

  if (width <= 0) {
    width = 2048;
  }
  if (height <= 0) {
    height = 2047;
  }
  if (iters <= 0) {
    iters = 256;
  }
  if (display <= 0) {
    display = 0;
  }

  unsigned *initial = (unsigned*)malloc(sizeof(unsigned) * height * width);

  // Initialize universe randomly
  for (unsigned y = 0; y < height; y++) {
    for (unsigned x = 0; x < width; x++) {
      initial[y * width + x] = (rand() % 3 == 0);
    }
  }

  double reference_time = 0;
  unsigned *reference;
  for (trials = 1; reference_time < TIMEOUT; trials *= 2) {
    // Warm-up
    reference = reference_life(height, width, initial, 1, display);
    free(reference);

    // Benchmark n runs of life
    reference_time = -wall_time();
    for (int i = 0; i < trials - 1; ++i){
      reference = reference_life(height, width, initial, iters, display);
      free(reference);
    }
    reference = reference_life(height, width, initial, iters, display);
    reference_time += wall_time();
  }
  trials /= 2;
  reference_time /= trials;

  double test_time = 0;
  unsigned *test;
  for (trials = 1; test_time < TIMEOUT; trials *= 2) {
    // Warm-up
    test = life(height, width, initial, 1, display);
    free(test);

    // Benchmark n runs of life
    test_time = -wall_time();
    for (int i = 0; i < trials - 1; ++i){
      test = life(height, width, initial, iters, display);
      free(test);
    }
    test = life(height, width, initial, iters, display);
    test_time += wall_time();
  }
  trials /= 2;
  test_time /= trials;

  // Ensure that reference and test results are identical
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

  printf("Reference time: %g\n", reference_time);
  printf("     Test time: %g\n", test_time);
  printf("Reference freq: %g\n", (height * width * iters) / reference_time);
  printf("     Test freq: %g\n", (height * width * iters) / test_time);
  printf("       Speedup: %g\n", reference_time / test_time);

  free(initial);
  free(reference);
  free(test);

  return 0;
}
