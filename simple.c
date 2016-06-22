#include <stdlib.h>

unsigned *life (const unsigned height,
                const unsigned width,
                const unsigned * const initial,
                const unsigned iters) {
  unsigned *universe = (unsigned*)malloc(sizeof(unsigned) * height * width);
  unsigned *new = (unsigned*)malloc(sizeof(unsigned) * height * width);
  for (unsigned y = 0; y < height; y++) {
    for (unsigned x = 0; x < width; x++) {
      universe[y * width + x] = initial[y * width + x];
    }
  }
  for (unsigned i = 0; i < iters; i ++) {
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        unsigned n = 0;
        for (int yy = y - 1; yy <= y + 1; yy++) {
          for (int xx = x - 1; xx <= x + 1; xx++) {
              //Directly add "universe" values to "n"
              n += universe[((yy + height) % height) * width
                            + ((xx + width) % width)];
            }
          }
          n -= universe[y * width + x];
          new[y * width + x] = (n == 3 || (n == 2 && universe[y * width + x]));
      }
    }
    //Instead of copying "new" into universe every time, just swap the pointers
    unsigned *tmp = universe;
    universe = new;
    new = tmp;
  }
  free(new);
  return universe;
}
