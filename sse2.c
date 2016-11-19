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
#include <stdint.h>
#include <pmmintrin.h>

#define             WORD (128/8)
#define        OUT_GHOST 1
#define      X_OUT_GHOST (((OUT_GHOST - 1)/WORD + 1) * WORD)
#define      Y_OUT_GHOST OUT_GHOST
#define         IN_GHOST (OUT_GHOST + 1)
#define       X_IN_GHOST ((OUT_GHOST/WORD + 1) * WORD)
#define       Y_IN_GHOST IN_GHOST
#define X_IN_GHOST_WORDS (X_IN_GHOST/WORD)

void *aligned_malloc(int size) {
    char *mem = malloc(sizeof(void*) + size + WORD - 1);
    void **ptr = (void**)(((uintptr_t)(mem + sizeof(void*) + WORD - 1)) & ~((uintptr_t)(WORD - 1)));
    ptr[-1] = mem;
    return ptr;
}

void aligned_free(void *ptr) {
    free(((void**)ptr)[-1]);
}

unsigned *life (const unsigned height,
                const unsigned width,
                const unsigned * const initial,
                const unsigned iters) {
  const unsigned padded_height = height + 2 * Y_IN_GHOST;
  const unsigned padded_width = width + 2 * X_IN_GHOST;
  const unsigned width_words = width/WORD;
  const unsigned padded_width_words = padded_width/WORD;
  uint8_t *universe = (uint8_t*)aligned_malloc(padded_height * padded_width);
  uint8_t *new = (uint8_t*)aligned_malloc(padded_height * padded_width);

  for (unsigned y = Y_IN_GHOST; y < height + Y_IN_GHOST; y++) {
    for (unsigned x = X_IN_GHOST; x < width + X_IN_GHOST; x++) {
      universe[(y * padded_width) + x] = initial[(y - Y_IN_GHOST) * width + x - X_IN_GHOST];
    }
  }

  for (unsigned i = 0; i < iters; i+= IN_GHOST) {

    //Because we assume the width is a multiple of the size of a SSE register,
    //we can use aligned loads and stores.
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

    for (unsigned j = 0; j < IN_GHOST & j + i < iters; j++) {
      //Set up a vector of ones
      const __m128i ones = _mm_set_epi8(1, 1, 1, 1, 1, 1, 1, 1,
                                        1, 1, 1, 1, 1, 1, 1, 1);
      //Set up a vector of twos
      const __m128i twos = _mm_slli_epi32(ones, 1);
      //Set up a vector of threes
      const __m128i threes = _mm_or_si128(ones, twos);
      for (unsigned y = (Y_IN_GHOST - Y_OUT_GHOST); y < height + Y_IN_GHOST + Y_OUT_GHOST; y++) {
        for (unsigned x = (X_IN_GHOST - X_OUT_GHOST); x + WORD <= width + X_IN_GHOST + X_OUT_GHOST; x += WORD) {
          __m128i n;
          __m128i alive;
          uint8_t *u = universe + (y - 1) * padded_width + x - 1;
          //This is an unaligned load
          n = _mm_lddqu_si128((__m128i*)u);
          n = _mm_add_epi8(_mm_load_si128((__m128i*)(u + 1)), n);
          n = _mm_add_epi8(_mm_lddqu_si128((__m128i*)(u + 2)), n);
          u += padded_width;
          n = _mm_add_epi8(_mm_lddqu_si128((__m128i*)u), n);
          //This is an aligned load
          alive = _mm_load_si128((__m128i*)(u + 1));
          n = _mm_add_epi8(_mm_lddqu_si128((__m128i*)(u + 2)), n);
          u += padded_width;
          n = _mm_add_epi8(_mm_lddqu_si128((__m128i*)u), n);
          n = _mm_add_epi8(_mm_load_si128((__m128i*)(u + 1)), n);
          n = _mm_add_epi8(_mm_lddqu_si128((__m128i*)(u + 2)), n);
          //The operation we are performing here is the same, but it looks
          //very different when written in SIMD instructions
          _mm_store_si128((__m128i*)(new + y * padded_width + x),
            _mm_or_si128(
            //We need to and with the ones vector here because the result of
            //comparison is either 0xFF or 0, and we need 1 or 0.
            _mm_and_si128(ones, _mm_cmpeq_epi8(n, threes)),
            _mm_and_si128(alive, _mm_cmpeq_epi8(n, twos))));
        }
      }
      uint8_t *tmp = universe;
      universe = new;
      new = tmp;
    }
  }
  unsigned *out = (unsigned*)malloc(sizeof(unsigned) * height * width);
  for (unsigned y = Y_IN_GHOST; y < height + Y_IN_GHOST; y++) {
    for (unsigned x = X_IN_GHOST; x < width + X_IN_GHOST; x++) {
      out[(y - Y_IN_GHOST) * width + x - X_IN_GHOST] = universe[(y * padded_width) + x];
    }
  }

  aligned_free(new);
  aligned_free(universe);
  return out;
}
