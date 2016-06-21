#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>

#define             WORD sizeof(unsigned)
#define        OUT_GHOST WORD * 1
#define         IN_GHOST (OUT_GHOST + 1)
#define       X_IN_GHOST ((OUT_GHOST/WORD + 1) * WORD) //should be multiple of word size
#define       Y_IN_GHOST IN_GHOST
#define X_IN_GHOST_WORDS (X_IN_GHOST/WORD)

unsigned *life (const unsigned height,
                const unsigned width,
                const unsigned * const initial,
                const unsigned iters,
                const unsigned display) {
  const unsigned padded_height = height + 2 * Y_IN_GHOST;
  const unsigned padded_width = width + 2 * X_IN_GHOST;
  const unsigned width_words = width/WORD;
  const unsigned padded_width_words = padded_width/WORD;
  uint8_t *universe = (uint8_t*)malloc(padded_height * padded_width);
  uint8_t *new = (uint8_t*)malloc(padded_height * padded_width);

  //pack into padded working array
  for (unsigned y = Y_IN_GHOST; y < height + Y_IN_GHOST; y++) {
    for (unsigned x = X_IN_GHOST; x < width + X_IN_GHOST; x++) {
      universe[(y * padded_width) + x] = initial[(y - Y_IN_GHOST) * width + x - X_IN_GHOST];
    }
  }
  for (unsigned i = 0; i < iters; i++) {
    //copy the ghost cells once every IN_GHOST iterations
    if (i % IN_GHOST == 0) {
      unsigned *universe_words = (unsigned*)universe;
      for (unsigned y = 0; y < padded_height; y++) {
        if (y < Y_IN_GHOST) {
          for (unsigned x = 0; x < X_IN_GHOST_WORDS; x++) {
            universe_words[y * padded_width_words + x] = universe_words[(y + height) * padded_width_words + x + width_words];
          }
          for (unsigned x = X_IN_GHOST_WORDS; x < width_words + X_IN_GHOST_WORDS; x++) {
            universe_words[y * padded_width_words + x] = universe_words[(y + height) * padded_width_words + x];
          }
          for (unsigned x = width_words + X_IN_GHOST_WORDS ; x < padded_width_words; x++) {
            universe_words[y * padded_width_words + x] = universe_words[(y + height) * padded_width_words + x - width_words];
          }
        } else if (y < height + Y_IN_GHOST) {
          for (unsigned x = 0; x < X_IN_GHOST_WORDS; x++) {
            universe_words[y * padded_width_words + x] = universe_words[y * padded_width_words + x + width_words];
          }
          for (unsigned x = width_words + X_IN_GHOST_WORDS ; x < padded_width_words; x++) {
            universe_words[y * padded_width_words + x] = universe_words[y * padded_width_words + x - width_words];
          }
        } else {
          for (unsigned x = 0; x < X_IN_GHOST_WORDS; x++) {
            universe_words[y * padded_width_words + x] = universe_words[(y - height) * padded_width_words + x + width_words];
          }
          for (unsigned x = X_IN_GHOST_WORDS; x < width_words + X_IN_GHOST_WORDS; x++) {
            universe_words[y * padded_width_words + x] = universe_words[(y - height) * padded_width_words + x];
          }
          for (unsigned x = width_words + X_IN_GHOST_WORDS ; x < padded_width_words; x++) {
            universe_words[y * padded_width_words + x] = universe_words[(y - height) * padded_width_words + x - width_words];
          }
        }
      }
    }
    //evolve
    for (unsigned y = (Y_IN_GHOST - OUT_GHOST); y < height + Y_IN_GHOST + OUT_GHOST; y++) {
      for (unsigned x = (X_IN_GHOST - OUT_GHOST); x < width + X_IN_GHOST + OUT_GHOST; x++) {
        unsigned n = 0;
        uint8_t *u = universe + (y - 1) * padded_width + x - 1;
        n += u[0];
        n += u[1];
        n += u[2];
        u += padded_width;
        n += u[0];
        unsigned alive = u[1];
        n += u[2];
        u += padded_width;
        n += u[0];
        n += u[1];
        n += u[2];
        new[y * padded_width + x] = (n == 3 || (n == 2 && alive));
      }
    }
    uint8_t *tmp = universe;
    universe = new;
    new = tmp;
  }
  //unpack into output array
  unsigned *out = (unsigned*)malloc(sizeof(unsigned) * height * width);
  for (unsigned y = Y_IN_GHOST; y < height + Y_IN_GHOST; y++) {
    for (unsigned x = X_IN_GHOST; x < width + X_IN_GHOST; x++) {
      out[(y - Y_IN_GHOST) * width + x - X_IN_GHOST] = universe[(y * padded_width) + x];
    }
  }
  free(new);
  free(universe);
  return out;
}
