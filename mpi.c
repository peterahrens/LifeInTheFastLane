#include <stdlib.h>
#include <stdint.h>
#include <immintrin.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

#define             WORD (256/8)
#define        OUT_GHOST 3
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
  int rank;
  int size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  //We will be arranging our processes in a grid. We are assuming that the
  //problem width is divisible by the square root of the number of processors
  //times the width of a word, and that the height is divisible by the
  //square root of the number of processes
  const unsigned side = (int)sqrt((double)size);
  const unsigned my_height = height/side;
  const unsigned my_width = width/side;
  const unsigned my_y = (rank / side);
  const unsigned my_x = (rank % side);

  const unsigned my_padded_height = my_height + 2 * Y_IN_GHOST;
  const unsigned my_padded_width = my_width + 2 * X_IN_GHOST;
  const unsigned my_width_words = my_width/WORD;
  const unsigned my_padded_width_words = my_padded_width/WORD;
  uint8_t *universe = (uint8_t*)aligned_malloc(my_padded_height * my_padded_width);
  uint8_t *new = (uint8_t*)aligned_malloc(my_padded_height * my_padded_width);
  //Moving these to the top makes the code faster. OpenMP Mysteries!
  const __m256i ones = _mm256_set_epi8(1, 1, 1, 1, 1, 1, 1, 1,
                                       1, 1, 1, 1, 1, 1, 1, 1,
                                       1, 1, 1, 1, 1, 1, 1, 1,
                                       1, 1, 1, 1, 1, 1, 1, 1);
  const __m256i twos = _mm256_slli_epi32(ones, 1);
  const __m256i threes = _mm256_or_si256(ones, twos);

  //Send the data to all the processes
  uint8_t *scatter_buffer_send;
  uint8_t *scatter_buffer_recv = (uint8_t*)aligned_malloc(my_height * my_width);
  if (rank == 0) {
    scatter_buf = (uint8_t*)aligned_malloc(height * width);
    for (unsigned their_y = 0; their_y < side; their_y++) {
      for (unsigned their_x = 0; their_x < side; their_x++) {
        for (unsigned y = 0; y < my_height; y++) {
          for (unsigned x = 0; x < my_width; x++) {
            scatter_buffer_send[(their_y * side + their_width) * (my_width + my_height) + y * my_width + x] =
              initial[(their_y * my_height + y) * width + their_x * my_width + x];
          }
        }
      }
    }
  }
  MPI_Scatter((const void*)scatter_buffer_send, my_height * my_width, MPI_UNSIGNED_CHAR,
              (const void*)scatter_buffer_recv, my_height * my_width, 0,
              MPI_COMM_WORLD);
  for (unsigned y = Y_IN_GHOST; y < Y_IN_GHOST + my_height; y++) {
    for (unsigned x = X_IN_GHOST; x < X_IN_GHOST + my_width; x++) {

      universe[(y * my_padded_width) + x] = scatter_buffer_recv[(y - Y_IN_GHOST) * my_width + x - X_IN_GHOST];
    }
  }

  unsigned     *ghost_buffer_top_left_send = aligned_malloc(X_IN_GHOST * Y_IN_GHOST);
  unsigned          *ghost_buffer_top_send = aligned_malloc(  my_width * Y_IN_GHOST);
  unsigned    *ghost_buffer_top_right_send = aligned_malloc(X_IN_GHOST * Y_IN_GHOST);
  unsigned        *ghost_buffer_right_send = aligned_malloc(X_IN_GHOST * my_height );
  unsigned *ghost_buffer_bottom_right_send = aligned_malloc(X_IN_GHOST * Y_IN_GHOST);
  unsigned       *ghost_buffer_bottom_send = aligned_malloc(  my_width * Y_IN_GHOST);
  unsigned  *ghost_buffer_bottom_left_send = aligned_malloc(X_IN_GHOST * Y_IN_GHOST);
  unsigned         *ghost_buffer_left_send = aligned_malloc(X_IN_GHOST * my_height );

  unsigned *ghost_buffer_bottom_right_recv = aligned_malloc(X_IN_GHOST * Y_IN_GHOST);
  unsigned       *ghost_buffer_bottom_recv = aligned_malloc(  my_width * Y_IN_GHOST);
  unsigned  *ghost_buffer_bottom_left_recv = aligned_malloc(X_IN_GHOST * Y_IN_GHOST);
  unsigned         *ghost_buffer_left_recv = aligned_malloc(X_IN_GHOST * my_height );
  unsigned     *ghost_buffer_top_left_recv = aligned_malloc(X_IN_GHOST * Y_IN_GHOST);
  unsigned          *ghost_buffer_top_recv = aligned_malloc(  my_width * Y_IN_GHOST);
  unsigned    *ghost_buffer_top_right_recv = aligned_malloc(X_IN_GHOST * Y_IN_GHOST);
  unsigned        *ghost_buffer_right_recv = aligned_malloc(X_IN_GHOST * my_height );

  MPI_Request top_left_req;
  MPI_Request top_req;
  MPI_Request top_right_req;
  MPI_Request right_req;
  MPI_Request bottom_right_req;
  MPI_Request bottom_req;
  MPI_Request bottom_left_req;
  MPI_Request left_req;

  for (unsigned i = 0; i < iters; i+= IN_GHOST) {
    __m256i *universe_words = (__m256i*)universe;

    //Now we send ghost zones to all of our neighbors.
    //Top Left
    for (unsigned y = 0; y < Y_IN_GHOST; y++) {
      for (unsigned x = 0; x < X_IN_GHOST_WORDS; x++) {
        _mm256_store_si256(ghost_buffer_top_left_send + y * X_IN_GHOST_WORDS + x,
          _mm256_load_si256(universe_words + (y + my_height) * my_padded_width_words + x + my_width_words));
      }
    }
    MPI_ISend(ghost_buffer_top_left_send, X_IN_GHOST * Y_IN_GHOST, MPI_UNSIGNED_CHAR, 
              ((my_y + side - 1) % side) * side + ((my_x + side - 1) % side), 0, MPI_COMM_WORLD,
              &top_left_req);

    //Top
    for (unsigned y = 0; y < Y_IN_GHOST; y++) {
      for (unsigned x = X_IN_GHOST_WORDS; x < X_IN_GHOST_WORDS + my_width; x++) {
        _mm256_store_si256(ghost_buffer_top_send + y * my_width + x - X_IN_GHOST_WORDS,
          _mm256_load_si256(universe_words + (y + my_height) * my_padded_width_words + x + my_width_words));
      }
    }
    MPI_ISend(ghost_buffer_top_send, X_IN_GHOST * Y_IN_GHOST, MPI_UNSIGNED_CHAR, 
              ((my_y + side - 1) % side) * side + my_x, 0, MPI_COMM_WORLD,
              &top_req);

    //Top Right
    for (unsigned y = 0; y < Y_IN_GHOST; y++) {
      for (unsigned x = X_IN_GHOST_WORDS + my_width; x < my_padded_width; x++) {
        _mm256_store_si256(ghost_buffer_top_right_send + y * X_IN_GHOST_WORDS + x - X_IN_GHOST_WORDS - my_width,
          _mm256_load_si256(universe_words + (y + my_height) * my_padded_width_words + x + my_width_words));
      }
    }
    MPI_ISend(ghost_buffer_top_send, X_IN_GHOST * Y_IN_GHOST, MPI_UNSIGNED_CHAR, 
              ((my_y + side - 1) % side) * side + ((my_x + 1) % side), 0, MPI_COMM_WORLD,
              &top_right_req);

    //Right
    for (unsigned y = Y_IN_GHOST; y < Y_IN_GHOST + my_height; y++) {
        for (unsigned x = X_IN_GHOST_WORDS + my_width; x < my_padded_width; x++) {
          _mm256_store_si256(ghost_buffer_right_send + y * X_IN_GHOST_WORDS + x - X_IN_GHOST_WORDS - my_width,
            _mm256_load_si256(universe_words + y * my_padded_width_words + x + my_width_words));
        }
    MPI_ISend(ghost_buffer_right_send, X_IN_GHOST * my_height, MPI_UNSIGNED_CHAR, 
              my_y * side + ((my_x + 1) % side), 0, MPI_COMM_WORLD,
              &right_req);

    //Bottom Right
    for (unsigned y = Y_IN_GHOST + my_height; y < my_padded_height; y++) {
      for (unsigned x = X_IN_GHOST_WORDS + my_width; x < my_padded_width; x++) {
        _mm256_store_si256(ghost_buffer_bottom_right_send + (y - Y_IN_GHOST - my_height) * X_IN_GHOST_WORDS + x - X_IN_GHOST_WORDS - my_width,
          _mm256_load_si256(universe_words + (y + my_height) * my_padded_width_words + x + my_width_words));
      }
    }
    MPI_ISend(ghost_buffer_bottom_send, X_IN_GHOST * Y_IN_GHOST, MPI_UNSIGNED_CHAR, 
              ((my_y + 1) % side) * side + ((my_x + 1) % side), 0, MPI_COMM_WORLD,
              &bottom_right_req);

    //Bottom
    for (unsigned y = Y_IN_GHOST + my_height; y < my_padded_height; y++) {
      for (unsigned x = X_IN_GHOST_WORDS; x < X_IN_GHOST_WORDS + my_width; x++) {
        _mm256_store_si256(ghost_buffer_bottom_send + (y - Y_IN_GHOST - my_height) * my_width + x - X_IN_GHOST_WORDS,
          _mm256_load_si256(universe_words + (y + my_height) * my_padded_width_words + x + my_width_words));
      }
    }
    MPI_ISend(ghost_buffer_bottom_send, X_IN_GHOST * Y_IN_GHOST, MPI_UNSIGNED_CHAR, 
              ((my_y + 1) % side) * side + my_x, 0, MPI_COMM_WORLD,
              &bottom_req);

    //Bottom Left
    for (unsigned y = Y_IN_GHOST + my_height; y < my_padded_height; y++) {
      for (unsigned x = 0; x < X_IN_GHOST_WORDS; x++) {
        _mm256_store_si256(ghost_buffer_bottom_left_send + (y - Y_IN_GHOST - my_height) * X_IN_GHOST_WORDS + x,
          _mm256_load_si256(universe_words + (y + my_height) * my_padded_width_words + x + my_width_words));
      }
    }
    MPI_ISend(ghost_buffer_bottom_left_send, X_IN_GHOST * Y_IN_GHOST, MPI_UNSIGNED_CHAR, 
              ((my_y + 1) % side) * side + ((my_x + side - 1) % side), 0, MPI_COMM_WORLD,
              &bottom_left_req);

    //Left
    for (unsigned y = Y_IN_GHOST; y < Y_IN_GHOST + my_height; y++) {
        for (unsigned x = 0; x < X_IN_GHOST_WORDS; x++) {
          _mm256_store_si256(ghost_buffer_left_send + y * X_IN_GHOST_WORDS + x,
            _mm256_load_si256(universe_words + y * my_padded_width_words + x + my_width_words));
        }
    MPI_ISend(ghost_buffer_left_send, X_IN_GHOST * my_height, MPI_UNSIGNED_CHAR, 
              my_y * side + ((my_x + side - 1) % side), 0, MPI_COMM_WORLD,
              &left_req);

    //Now we recieve ghost zones from all of our neighbors.
    //Top Left
    MPI_Recv(ghost_buffer_top_left_recv, X_IN_GHOST * Y_IN_GHOST, MPI_UNSIGNED_CHAR, 
             ((my_y + side - 1) % side) * side + ((my_x + side - 1) % side), 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    for (unsigned y = 0; y < Y_IN_GHOST; y++) {
      for (unsigned x = 0; x < X_IN_GHOST_WORDS; x++) {
        _mm256_store_si256(universe_words + (y + my_height) * my_padded_width_words + x + my_width_words,
          _mm256_load_si256(ghost_buffer_top_left_recv + y * X_IN_GHOST_WORDS + x));
      }
    }

    //Top
    MPI_Recv(ghost_buffer_top_recv, X_IN_GHOST * Y_IN_GHOST, MPI_UNSIGNED_CHAR, 
             ((my_y + side - 1) % side) * side + my_x, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    for (unsigned y = 0; y < Y_IN_GHOST; y++) {
      for (unsigned x = X_IN_GHOST_WORDS; x < X_IN_GHOST_WORDS + my_width; x++) {
        _mm256_store_si256(universe_words + (y + my_height) * my_padded_width_words + x + my_width_words,
          _mm256_load_si256(ghost_buffer_top_recv + y * my_width + x - X_IN_GHOST_WORDS));
      }
    }

    //Top Right
    MPI_Recv(ghost_buffer_top_recv, X_IN_GHOST * Y_IN_GHOST, MPI_UNSIGNED_CHAR, 
             ((my_y + side - 1) % side) * side + ((my_x + 1) % side), 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    for (unsigned y = 0; y < Y_IN_GHOST; y++) {
      for (unsigned x = X_IN_GHOST_WORDS + my_width; x < my_padded_width; x++) {
        _mm256_store_si256(universe_words + (y + my_height) * my_padded_width_words + x + my_width_words,
          _mm256_load_si256(ghost_buffer_top_right_recv + y * X_IN_GHOST_WORDS + x - X_IN_GHOST_WORDS - my_width));
      }
    }

    //Right
    MPI_Recv(ghost_buffer_right_recv, X_IN_GHOST * my_height, MPI_UNSIGNED_CHAR, 
             my_y * side + ((my_x + 1) % side), 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    for (unsigned y = Y_IN_GHOST; y < Y_IN_GHOST + my_height; y++) {
        for (unsigned x = X_IN_GHOST_WORDS + my_width; x < my_padded_width; x++) {
          _mm256_store_si256(universe_words + y * my_padded_width_words + x + my_width_words,
            _mm256_load_si256(ghost_buffer_right_recv + y * X_IN_GHOST_WORDS + x - X_IN_GHOST_WORDS - my_width));
        }
    }

    //Bottom Right
    MPI_Recv(ghost_buffer_bottom_recv, X_IN_GHOST * Y_IN_GHOST, MPI_UNSIGNED_CHAR, 
             ((my_y + 1) % side) * side + ((my_x + 1) % side), 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    for (unsigned y = Y_IN_GHOST + my_height; y < my_padded_height; y++) {
      for (unsigned x = X_IN_GHOST_WORDS + my_width; x < my_padded_width; x++) {
        _mm256_store_si256(universe_words + (y + my_height) * my_padded_width_words + x + my_width_words,
          _mm256_load_si256(ghost_buffer_bottom_right_recv + (y - Y_IN_GHOST - my_height) * X_IN_GHOST_WORDS + x - X_IN_GHOST_WORDS - my_width));
      }
    }

    //Bottom
    MPI_Recv(ghost_buffer_bottom_recv, X_IN_GHOST * Y_IN_GHOST, MPI_UNSIGNED_CHAR, 
             ((my_y + 1) % side) * side + my_x, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    for (unsigned y = Y_IN_GHOST + my_height; y < my_padded_height; y++) {
      for (unsigned x = X_IN_GHOST_WORDS; x < X_IN_GHOST_WORDS + my_width; x++) {
        _mm256_store_si256(universe_words + (y + my_height) * my_padded_width_words + x + my_width_words,
          _mm256_load_si256(ghost_buffer_bottom_recv + (y - Y_IN_GHOST - my_height) * my_width + x - X_IN_GHOST_WORDS));
      }
    }

    //Bottom Left
    MPI_Recv(ghost_buffer_bottom_left_recv, X_IN_GHOST * Y_IN_GHOST, MPI_UNSIGNED_CHAR, 
             ((my_y + 1) % side) * side + ((my_x + side - 1) % side), 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    for (unsigned y = Y_IN_GHOST + my_height; y < my_padded_height; y++) {
      for (unsigned x = 0; x < X_IN_GHOST_WORDS; x++) {
        _mm256_store_si256(universe_words + (y + my_height) * my_padded_width_words + x + my_width_words,
          _mm256_load_si256(ghost_buffer_bottom_left_recv + (y - Y_IN_GHOST - my_height) * X_IN_GHOST_WORDS + x));
      }
    }

    //Left
    MPI_Recv(ghost_buffer_left_recv, X_IN_GHOST * my_height, MPI_UNSIGNED_CHAR, 
             my_y * side + ((my_x + side - 1) % side), 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    for (unsigned y = Y_IN_GHOST; y < Y_IN_GHOST + my_height; y++) {
        for (unsigned x = 0; x < X_IN_GHOST_WORDS; x++) {
          _mm256_store_si256(universe_words + y * my_padded_width_words + x + my_width_words,
            _mm256_load_si256(ghost_buffer_left_recv + y * X_IN_GHOST_WORDS + x));
        }
    }

    #pragma omp parallel
    {
      uint8_t *my_universe = universe;
      uint8_t *my_new = new;

      for (unsigned j = 0; j < IN_GHOST & j + i < iters; j++) {
        #pragma omp for
        for (unsigned y = (Y_IN_GHOST - Y_OUT_GHOST); y < my_height + Y_IN_GHOST + Y_OUT_GHOST; y++) {
          for (unsigned x = (X_IN_GHOST - X_OUT_GHOST); x + WORD <= my_width + X_IN_GHOST + X_OUT_GHOST; x += WORD) {
            __m256i n;
            __m256i alive;
            uint8_t *u = my_universe + (y - 1) * my_padded_width + x - 1;
            n = _mm256_lddqu_si256((__m256i*)u);
            n = _mm256_add_epi8(_mm256_load_si256((__m256i*)(u + 1)), n);
            n = _mm256_add_epi8(_mm256_lddqu_si256((__m256i*)(u + 2)), n);
            u += my_padded_width;
            n = _mm256_add_epi8(_mm256_lddqu_si256((__m256i*)u), n);
            alive = _mm256_load_si256((__m256i*)(u + 1));
            n = _mm256_add_epi8(_mm256_lddqu_si256((__m256i*)(u + 2)), n);
            u += my_padded_width;
            n = _mm256_add_epi8(_mm256_lddqu_si256((__m256i*)u), n);
            n = _mm256_add_epi8(_mm256_load_si256((__m256i*)(u + 1)), n);
            n = _mm256_add_epi8(_mm256_lddqu_si256((__m256i*)(u + 2)), n);
            _mm256_store_si256((__m256i*)(my_new + y * my_padded_width + x),
              _mm256_or_si256(
              _mm256_and_si256(ones, _mm256_cmpeq_epi8(n, threes)),
              _mm256_and_si256(alive, _mm256_cmpeq_epi8(n, twos))));
          }
        }
        uint8_t *tmp = my_universe;
        my_universe = my_new;
        my_new = tmp;
      }
      #pragma omp single
      {
        universe = my_universe;
        new = my_new;
      }
    }

    //Before we start another iteration and start sending again, let's make sure
    //that everyone has recieved our messages.
    MPI_Wait(&top_left_req, MPI_STATUS_IGNORE);
    MPI_Wait(&top_req, MPI_STATUS_IGNORE);
    MPI_Wait(&top_right_req, MPI_STATUS_IGNORE);
    MPI_Wait(&right_req, MPI_STATUS_IGNORE);
    MPI_Wait(&bottom_right_req, MPI_STATUS_IGNORE);
    MPI_Wait(&bottom_req, MPI_STATUS_IGNORE);
    MPI_Wait(&bottom_left_req, MPI_STATUS_IGNORE);
    MPI_Wait(&left_req, MPI_STATUS_IGNORE);
  }

  unsigned *out = NULL;
  //Recieve the data from all the processes
  for (unsigned y = Y_IN_GHOST; y < Y_IN_GHOST + my_height; y++) {
    for (unsigned x = X_IN_GHOST; x < X_IN_GHOST + my_width; x++) {
      scatter_buffer_recv[(y - Y_IN_GHOST) * my_width + x - X_IN_GHOST] = universe[(y * my_padded_width) + x];
    }
  }
  MPI_Gather((const void*)scatter_buffer_recv, my_height * my_width, MPI_UNSIGNED_CHAR,
             (const void*)scatter_buffer_send, my_height * my_width, 0,
             MPI_COMM_WORLD);
  if (rank == 0) {
    out = (unsigned*)malloc(sizeof(unsigned) * my_height * my_width);
    for (unsigned their_y = 0; their_y < side; their_y++) {
      for (unsigned their_x = 0; their_x < side; their_x++) {
        for (unsigned y = 0; y < my_height; y++) {
          for (unsigned x = 0; x < my_width; x++) {
            out[(their_y * my_height + y) * width + their_x * my_width + x] =
            scatter_buffer_send[(their_y * side + their_width) * (my_width + my_height) + y * my_width + x];
          }
        }
      }
    }
  }

  aligned_free(new);
  aligned_free(universe);
  aligned_free(ghost_buffer_top_left_send);
  aligned_free(ghost_buffer_top_send);
  aligned_free(ghost_buffer_top_right_send);
  aligned_free(ghost_buffer_right_send);
  aligned_free(ghost_buffer_bottom_right_send);
  aligned_free(ghost_buffer_bottom_send);
  aligned_free(ghost_buffer_bottom_left_send);
  aligned_free(ghost_buffer_left_send);
  aligned_free(ghost_buffer_top_left_recv);
  aligned_free(ghost_buffer_top_recv);
  aligned_free(ghost_buffer_top_right_recv);
  aligned_free(ghost_buffer_right_recv);
  aligned_free(ghost_buffer_bottom_right_recv);
  aligned_free(ghost_buffer_bottom_recv);
  aligned_free(ghost_buffer_bottom_left_recv);
  aligned_free(ghost_buffer_left_recv);
  if (rank == 0) {
    aligned_free(scatter_buffer_send);
  }
  aligned_free(scatter_buffer_recv);
  return out;
}
