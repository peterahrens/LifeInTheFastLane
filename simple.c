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
