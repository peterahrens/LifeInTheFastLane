#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <pmmintrin.h>

#define             WORD (128/8)
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
    if (i % IN_GHOST == 0){
      __m128i *universe_words = (__m128i*)universe;
      for (unsigned y = 0; y < padded_height; y++) {
        if (y < Y_IN_GHOST) {
          for (unsigned x = 0; x < X_IN_GHOST_WORDS; x++) {
            _mm_store_si128(universe_words + y * padded_width_words + x,
              _mm_load_si128(universe_words + (y + height) * padded_width_words + x + width_words));
          }
          for (unsigned x = X_IN_GHOST_WORDS; x < width_words + X_IN_GHOST_WORDS; x++) {
            _mm_store_si128(universe_words + y * padded_width_words + x,
              _mm_load_si128(universe_words + (y + height) * padded_width_words + x));
          }
          for (unsigned x = width_words + X_IN_GHOST_WORDS ; x < padded_width_words; x++) {
            _mm_store_si128(universe_words + y * padded_width_words + x,
              _mm_load_si128(universe_words + (y + height) * padded_width_words + x - width_words));
          }
        } else if (y < height + Y_IN_GHOST) {
          for (unsigned x = 0; x < X_IN_GHOST_WORDS; x++) {
            _mm_store_si128(universe_words + y * padded_width_words + x,
              _mm_load_si128(universe_words + y * padded_width_words + x + width_words));
          }
          for (unsigned x = width_words + X_IN_GHOST_WORDS ; x < padded_width_words; x++) {
            _mm_store_si128(universe_words + y * padded_width_words + x,
              _mm_load_si128(universe_words + y * padded_width_words + x - width_words));
          }
        } else {
          for (unsigned x = 0; x < X_IN_GHOST_WORDS; x++) {
            _mm_store_si128(universe_words + y * padded_width_words + x,
              _mm_load_si128(universe_words + (y - height) * padded_width_words + x + width_words));
          }
          for (unsigned x = X_IN_GHOST_WORDS; x < width_words + X_IN_GHOST_WORDS; x++) {
            _mm_store_si128(universe_words + y * padded_width_words + x,
              _mm_load_si128(universe_words + (y - height) * padded_width_words + x));
          }
          for (unsigned x = width_words + X_IN_GHOST_WORDS ; x < padded_width_words; x++) {
            _mm_store_si128(universe_words + y * padded_width_words + x,
              _mm_load_si128(universe_words + (y - height) * padded_width_words + x - width_words));
          }
        }
      }
    }
    //evolve
    __m128i ones = _mm_set_epi8(1, 1, 1, 1, 1, 1, 1, 1,
                                1, 1, 1, 1, 1, 1, 1, 1);
    __m128i twos = _mm_slli_epi32(ones, 1);
    __m128i threes = _mm_or_si128(ones, twos);
    for (unsigned y = (Y_IN_GHOST - OUT_GHOST); y < height + Y_IN_GHOST + OUT_GHOST; y++) {
      for (unsigned x = (X_IN_GHOST - OUT_GHOST); x + WORD <= width + X_IN_GHOST + OUT_GHOST; x += WORD) {
        __m128i n;
        __m128i alive;
        uint8_t *u = universe + (y - 1) * padded_width + x - 1;
        n = _mm_loadu_si128((__m128i*)u);
        n = _mm_add_epi8(_mm_load_si128((__m128i*)(u + 1)), n);
        n = _mm_add_epi8(_mm_loadu_si128((__m128i*)(u + 2)), n);
        u += padded_width;
        n = _mm_add_epi8(_mm_loadu_si128((__m128i*)u), n);
        alive = _mm_load_si128((__m128i*)(u + 1));
        n = _mm_add_epi8(_mm_loadu_si128((__m128i*)(u + 2)), n);
        u += padded_width;
        n = _mm_add_epi8(_mm_loadu_si128((__m128i*)u), n);
        n = _mm_add_epi8(_mm_load_si128((__m128i*)(u + 1)), n);
        n = _mm_add_epi8(_mm_loadu_si128((__m128i*)(u + 2)), n);
        _mm_store_si128((__m128i*)(new + y * padded_width + x),
          _mm_or_si128(
          _mm_and_si128(ones, _mm_cmpeq_epi8(n, threes)),
          _mm_and_si128(alive, _mm_cmpeq_epi8(n, twos))));
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
