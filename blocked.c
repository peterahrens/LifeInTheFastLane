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

#define             WORD sizeof(unsigned)
#define        OUT_GHOST 0
#define         IN_GHOST (OUT_GHOST + 1)
#define       X_IN_GHOST ((OUT_GHOST/WORD + 1) * WORD)
#define       Y_IN_GHOST IN_GHOST
#define X_IN_GHOST_WORDS (X_IN_GHOST/WORD)
#define          X_BLOCK WORD * 256
#define          Y_BLOCK 256

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

  for (unsigned i = 0; i < iters; i += IN_GHOST) {

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

    for (unsigned j = 0; j < IN_GHOST && i + j < iters; j++) {
      //Now the outer loops progress block by block.
      for (unsigned y = (Y_IN_GHOST - OUT_GHOST); y < height + Y_IN_GHOST + OUT_GHOST; y += Y_BLOCK) {
        for (unsigned x = (X_IN_GHOST - OUT_GHOST); x < width + X_IN_GHOST + OUT_GHOST; x += X_BLOCK) {
          //The inner loops progress one by one.
          for (unsigned yy = y; yy < y + Y_BLOCK && yy < height + Y_IN_GHOST + OUT_GHOST; yy++) {
            for (unsigned xx = x; xx < x + X_BLOCK && xx < width + X_IN_GHOST + OUT_GHOST; xx++) {
              unsigned n = 0;
              uint8_t *u = universe + (yy - 1) * padded_width + xx - 1;
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
              new[yy * padded_width + xx] = (n == 3 || (n == 2 && alive));
            }
          }
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
