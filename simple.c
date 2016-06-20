#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

static void evolve_into (const unsigned height,
                         const unsigned width,
                         unsigned * const universe,
                         unsigned * const new) {
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
}

unsigned *life (const unsigned height,
                const unsigned width,
                const unsigned * const initial,
                const unsigned iters,
                const unsigned display) {
  unsigned *universe = (unsigned*)malloc(sizeof(unsigned) * height * width);
  unsigned *new = (unsigned*)malloc(sizeof(unsigned) * height * width);
  for (unsigned y = 0; y < height; y++) {
    for (unsigned x = 0; x < width; x++) {
      universe[y * width + x] = initial[y * width + x];
    }
  }
  for (unsigned i = 0; i < iters; i ++) {
    evolve_into(height, width, universe, new);
    unsigned *tmp = universe;
    universe = new;
    new = tmp;
  }
  free(new);
  return universe;
}
