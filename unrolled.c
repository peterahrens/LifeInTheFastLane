#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <immintrin.h>

#define             WORD (256/8)
#define        OUT_GHOST WORD * 1
#define         IN_GHOST (OUT_GHOST + 1)
#define       X_IN_GHOST ((OUT_GHOST/WORD + 1) * WORD) //should be multiple of word size
#define       Y_IN_GHOST IN_GHOST
#define X_IN_GHOST_WORDS (X_IN_GHOST/WORD)

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

static void show_char (unsigned height, unsigned width, uint8_t *universe) {
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
      __m256i *universe_words = (__m256i*)universe;
      for (unsigned y = 0; y < padded_height; y++) {
        if (y < Y_IN_GHOST) {
          for (unsigned x = 0; x < X_IN_GHOST_WORDS; x++) {
            _mm256_store_si256(universe_words + y * padded_width_words + x,
              _mm256_load_si256(universe_words + (y + height) * padded_width_words + x + width_words));
          }
          for (unsigned x = X_IN_GHOST_WORDS; x < width_words + X_IN_GHOST_WORDS; x++) {
            _mm256_store_si256(universe_words + y * padded_width_words + x,
              _mm256_load_si256(universe_words + (y + height) * padded_width_words + x));
          }
          for (unsigned x = width_words + X_IN_GHOST_WORDS ; x < padded_width_words; x++) {
            _mm256_store_si256(universe_words + y * padded_width_words + x,
              _mm256_load_si256(universe_words + (y + height) * padded_width_words + x - width_words));
          }
        } else if (y < height + Y_IN_GHOST) {
          for (unsigned x = 0; x < X_IN_GHOST_WORDS; x++) {
            _mm256_store_si256(universe_words + y * padded_width_words + x,
              _mm256_load_si256(universe_words + y * padded_width_words + x + width_words));
          }
          for (unsigned x = width_words + X_IN_GHOST_WORDS ; x < padded_width_words; x++) {
            _mm256_store_si256(universe_words + y * padded_width_words + x,
              _mm256_load_si256(universe_words + y * padded_width_words + x - width_words));
          }
        } else {
          for (unsigned x = 0; x < X_IN_GHOST_WORDS; x++) {
            _mm256_store_si256(universe_words + y * padded_width_words + x,
              _mm256_load_si256(universe_words + (y - height) * padded_width_words + x + width_words));
          }
          for (unsigned x = X_IN_GHOST_WORDS; x < width_words + X_IN_GHOST_WORDS; x++) {
            _mm256_store_si256(universe_words + y * padded_width_words + x,
              _mm256_load_si256(universe_words + (y - height) * padded_width_words + x));
          }
          for (unsigned x = width_words + X_IN_GHOST_WORDS ; x < padded_width_words; x++) {
            _mm256_store_si256(universe_words + y * padded_width_words + x,
              _mm256_load_si256(universe_words + (y - height) * padded_width_words + x - width_words));
          }
        }
      }
    }
    //evolve
    __m256i ones = _mm256_set_epi8(1, 1, 1, 1, 1, 1, 1, 1,
                                   1, 1, 1, 1, 1, 1, 1, 1,
                                   1, 1, 1, 1, 1, 1, 1, 1,
                                   1, 1, 1, 1, 1, 1, 1, 1);
    __m256i twos = _mm256_slli_epi32(ones, 1);
    __m256i threes = _mm256_or_si256(ones, twos);
    for (unsigned y = (Y_IN_GHOST - OUT_GHOST); y + 2 <= height + Y_IN_GHOST + OUT_GHOST; y += 2) {
      for (unsigned x = (X_IN_GHOST - OUT_GHOST); x + WORD <= width + X_IN_GHOST + OUT_GHOST; x += WORD) {
        __m256i n_0;
        __m256i t;
        __m256i n_1;
        __m256i alive_0;
        __m256i alive_1;
        uint8_t *u = universe + (y - 1) * padded_width + x - 1;
        //n_0 gets (0, 0), (0, 1), (0, 2)
        n_0 = _mm256_loadu_si256((__m256i*)u);
        n_0 = _mm256_add_epi8(_mm256_load_si256((__m256i*)(u + 1)), n_0);
        n_0 = _mm256_add_epi8(_mm256_loadu_si256((__m256i*)(u + 2)), n_0);
        u += padded_width;
        //n_0 gets (1, 0), (1, 2)
        //n_1 gets (1, 0), (1, 1), (1, 2)
        n_1 = _mm256_loadu_si256((__m256i*)u);
        alive_0 = _mm256_load_si256((__m256i*)(u + 1));
        n_1 = _mm256_add_epi8(_mm256_loadu_si256((__m256i*)(u + 2)), n_1);
        n_0 = _mm256_add_epi8(n_1, n_0);
        n_1 = _mm256_add_epi8(alive_0, n_1);
        u += padded_width;
        //n_0 gets (2, 0), (2, 1), (2, 2)
        //n_1 gets (2, 0), (2, 2)
        t = _mm256_loadu_si256((__m256i*)u);
        alive_1 = _mm256_load_si256((__m256i*)(u + 1));
        t = _mm256_add_epi8(_mm256_loadu_si256((__m256i*)(u + 2)), t);
        n_0 = _mm256_add_epi8(t, n_0);
        n_0 = _mm256_add_epi8(alive_1, n_0);
        n_1 = _mm256_add_epi8(t, n_1);
        u += padded_width;
        //n_1 gets (3, 0), (3, 1), (3, 2)
        n_1 = _mm256_add_epi8(_mm256_loadu_si256((__m256i*)u), n_1);
        n_1 = _mm256_add_epi8(_mm256_load_si256((__m256i*)(u + 1)), n_1);
        n_1 = _mm256_add_epi8(_mm256_loadu_si256((__m256i*)(u + 2)), n_1);

        _mm256_store_si256((__m256i*)(new + y * padded_width + x),
          _mm256_or_si256(
          _mm256_and_si256(ones, _mm256_cmpeq_epi8(n_0, threes)),
          _mm256_and_si256(alive_0, _mm256_cmpeq_epi8(n_0, twos))));
        _mm256_store_si256((__m256i*)(new + (y + 1) * padded_width + x),
          _mm256_or_si256(
          _mm256_and_si256(ones, _mm256_cmpeq_epi8(n_1, threes)),
          _mm256_and_si256(alive_1, _mm256_cmpeq_epi8(n_1, twos))));
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
