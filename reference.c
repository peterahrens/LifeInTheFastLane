#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>

static void show (unsigned height, unsigned width, unsigned *universe) {
  printf("\033[H");
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
        printf(universe[y * width + x] ? "\033[07m  \033[m" : "  ");
    }
    printf("\033[E");
  }
  fflush(stdout);
}

static void evolve (unsigned height,
                    unsigned width,
                    unsigned *universe,
                    unsigned *new) {
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      unsigned n = 0;
      for (int yy = y - 1; yy <= y + 1; yy++) {
        for (int xx = x - 1; xx <= x + 1; xx++) {
          if (universe[((yy + height) % height) * width +
              ((xx + width) % width)]) {
            n++;
          }
        }
      }
      if (universe[y * width + x]) {
        n--;
      }
      new[y * width + x] = (n == 3 || (n == 2 && universe[y * width + x]));
    }
  }
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      universe[y * width + x] = new[y * width + x];
    }
  }
}

unsigned *reference_life (unsigned height,
                          unsigned width,
                          unsigned *initial,
                          unsigned iters,
                          unsigned display) {
  unsigned *universe = (unsigned*)malloc(sizeof(unsigned) * height * width);
  unsigned *new = (unsigned*)malloc(sizeof(unsigned) * height * width);
  for (unsigned y = 0; y < height; y++) {
    for (unsigned x = 0; x < width; x++) {
      universe[y * width + x] = initial[y * width + x];
    }
  }
  for (unsigned i = 0; i < iters; i++) {
    if (display && i % display == 0) {
      show(height, width, universe);
    }
    evolve(height, width, universe, new);
  }
  if (display) {
    show(height, width, universe);
  }
  free(new);
  return universe;
}
