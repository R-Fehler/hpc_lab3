#include <mpi.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// OPTIONAL: comment this out for console output
//#define CONSOLE_OUTPUT

#define calcIndex(width, x, y) ((y) * (width) + (x))
#define ALIVE 1
#define DEAD 0

#define X 0  //
#define Y 1  //

int rank = 0;       // The current MPI rank in the global communicator.
int rank_cart = 0;  // The current MPI rank in the cart communicator.
int num_tasks;      // The number of processes
int coords[2];

// evolve
int starts[2];
int ends[2];

MPI_Datatype filetype;  //
MPI_Comm cart_comm;     // Communicator for the cartesian grid
MPI_File file;          // A shared file pointer
MPI_Datatype memtype;   // A new created type for the inner and outer data
                        // including the ghost layer

void myexit(const char *s, ...) {
  va_list args;
  va_start(args, s);
  if (rank == 0) {
    vprintf(s, args);
    printf("\n");
  }
  va_end(args);
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

char vtk_header[2048];
void create_vtk_header(char *header, int width, int height, int timestep) {
  char buffer[1024];
  header[0] = '\0';
  strcat(header, "# vtk DataFile Version 3.0\n");
  snprintf(buffer, sizeof(buffer), "Gameoflife timestep %d \n", timestep);
  strcat(header, buffer);
  strcat(header, "BINARY\n");
  strcat(header, "DATASET STRUCTURED_POINTS\n");
  snprintf(buffer, sizeof(buffer), "DIMENSIONS %d %d 1\n", width, height);
  strcat(header, buffer);
  strcat(header, "SPACING 1.0 1.0 1.0\n");
  strcat(header, "ORIGIN 0 0 0\n");
  snprintf(buffer, sizeof(buffer), "POINT_DATA %ld\n", width * height);
  strcat(header, buffer);
  strcat(header, "SCALARS data char 1\n");
  strcat(header, "LOOKUP_TABLE default\n");
}

void write_field(char *currentfield, int width, int height, int timestep) {
  if (timestep == 0) {
    if (rank_cart == 0) {
      mkdir("./gol/", 0777);
    }
    create_vtk_header(vtk_header, width, height, timestep);
  }

  printf("writing timestep %d\n", timestep);
  char filename[1024];
  snprintf(filename, 1024, "./gol/gol-%05d.vtk", timestep);
  MPI_Offset header_offset = (MPI_Offset)strlen(vtk_header);

  /* TODO Create a new file handle for collective I/O
   *      Use the global 'file' variable.
   */
  MPI_File_open(cart_comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &file);
  if (rank_cart == 0) {
    MPI_File_write(file, vtk_header, strlen(vtk_header), MPI_CHAR,
                   MPI_STATUS_IGNORE);
  }
  /* TODO Set the file view for the file handle using collective I/O
   *
   */
  MPI_File_set_view(file, header_offset, MPI_CHAR, filetype, "native",
                    MPI_INFO_NULL);
  // rc = ...

  /* TODO Write the data using collective I/O
   *
   */
  MPI_File_write_all(file, currentfield, width * height, MPI_CHAR,
                     MPI_STATUS_IGNORE);

  /* TODO Close the file handle.
   *
   */
  MPI_File_close(&file);
}

void evolve(char *currentfield, char *newfield, int width, int height) {
  // TODO traverse through each voxel and implement game of live logic
  for (int y = starts[Y]; y < ends[Y]; y++) {
    int yi = y * width;
    char *ci = currentfield + yi;
    char *ni = newfield + yi;

    for (int x = starts[X]; x < ends[X]; x++) {
      // printf("%d;%d;%d\n" , omp_get_thread_num(),x, y);
      int neighbours = 0;

      for (int ni = 0; ni < 9; ni++) {
        if (ni == 4) {
          continue;
        }
        int ny = (ni / 3) - 1;
        int nx = (ni % 3) - 1;

        char *nyi = ci + x + (ny * width);
        char neighbourCell = *(nyi + nx);
        neighbours += neighbourCell;
      }

      char cell = ci[x];
      switch (cell) {
        case ALIVE:
          if (neighbours < 2) {
            ni[x] = DEAD;
          } else if (neighbours > 3) {
            ni[x] = DEAD;
          } else {
            ni[x] = ALIVE;
          }
          break;
        case DEAD:
          if (neighbours == 3) {
            ni[x] = ALIVE;
          } else {
            ni[x] = DEAD;
          }
      }
    }
  }
  // HINT: avoid boundaries
}

void filling_random(char *currentfield, int width, int height) {
  int i;
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      i = calcIndex(width, x, y);
      currentfield[i] =
          (rand() < RAND_MAX / 10) ? 1 : 0;  ///< init domain randomly
    }
  }
}

void filling_runner(char *currentfield, int width, int height) {
  currentfield[calcIndex(width, width / 2 + 0, height / 2 + 1)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 1, height / 2 + 2)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 2, height / 2 + 0)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 2, height / 2 + 1)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 2, height / 2 + 2)] = ALIVE;
}

void exchange(char *currentfield, int width, int height) {
  int x = coords[X];
  int y = coords[Y];
  int n[2];
  n[X] = x;
  n[Y] = y - 1;
  int s[2];
  s[X] = x;
  s[Y] = y + 1;
  int w[2];
  w[X] = x - 1;
  w[Y] = y;
  int e[2];
  e[X] = x + 1;
  e[Y] = y;

  // NORTH

  int n_id;
  MPI_Cart_rank(cart_comm, n, &n_id);
  char *n_sendbuf = currentfield + width;
  char *n_recvbuf = currentfield;
  MPI_Request n_ignore;

  // delete
  // MPI_ISend(n_sendbuf, width, MPI_CHAR, n_id, 0, cart_comm, &n_ignore);

  // SOUTH

  int s_id;
  MPI_Cart_rank(cart_comm, s, &s_id);
  char *s_sendbuf = currentfield + width * (height - 2);
  char *s_recvbuf = currentfield + width * (height - 1);
  MPI_Request s_ignore;
  // MPI_ISend(s_sendbuf, width, MPI_CHAR, s_id, 1, cart_comm, &s_ignore);

  // MPI_Recv(n_recvbuf, width, MPI_CHAR, n_id, 1, cart_comm,
  // MPI_STATUS_IGNORE); MPI_Recv(s_recvbuf, width, MPI_CHAR, s_id, 0,
  // cart_comm, MPI_STATUS_IGNORE);

  // WEST
  char w_sendbuf[height];

  int w_id;
  MPI_Cart_rank(cart_comm, w, &w_id);
  char w_recvbuf[height];
  for (int y = 0; y < height; y++) {
    w_recvbuf[y] = 0;
  }

  // EAST
  char e_sendbuf[height];

  int e_id;
  MPI_Cart_rank(cart_comm, e, &e_id);
  char e_recvbuf[height];
  for (int y = 0; y < height; y++) {
    e_recvbuf[y] = 0;
  }

  // SEND

  printf("Exchanging north bounds %d\n", rank_cart);
  printf("fill north %d ", n_sendbuf);
  for (int y = 0; y < width; y++) {
    printf("%d ", n_sendbuf[y]);
  }
  printf("\nsend north %d to %d\n", n_sendbuf, s_recvbuf);
  MPI_Sendrecv(n_sendbuf, width, MPI_CHAR, n_id, 0, s_recvbuf, width, MPI_CHAR,
               s_id, 0, cart_comm, MPI_STATUS_IGNORE);
  printf("recv south %d ", s_recvbuf);
  for (int y = 0; y < width; y++) {
    printf("%d ", s_recvbuf[y]);
  }
  printf("\n");

  printf("Exchanging south bounds %d\n", rank_cart);
  printf("fill south %d ", s_sendbuf);
  for (int y = 0; y < width; y++) {
    printf("%d ", s_sendbuf[y]);
  }
  printf("\nsend south %d to %d\n", s_sendbuf, n_recvbuf);
  MPI_Sendrecv(s_sendbuf, width, MPI_CHAR, s_id, 1, n_recvbuf, width, MPI_CHAR,
               n_id, 1, cart_comm, MPI_STATUS_IGNORE);
  printf("recv north %d ", n_recvbuf);
  for (int y = 0; y < width; y++) {
    printf("%d ", n_recvbuf[y]);
  }
  printf("\n");

  printf("Exchanging west bounds %d\n", rank_cart);
  printf("fill west %d ", w_recvbuf);
  for (int y = 0; y < height; y++) {
    w_sendbuf[y] = currentfield[calcIndex(width, 1, y)];
    printf("%d ", w_sendbuf[y]);
  }
  printf("\nsend west %d to %d\n", w_sendbuf, e_recvbuf);
  MPI_Sendrecv(w_sendbuf, width, MPI_CHAR, w_id, 3, e_recvbuf, width, MPI_CHAR,
               e_id, 3, cart_comm, MPI_STATUS_IGNORE);
  printf("recv east %d ", e_recvbuf);
  for (int y = 0; y < height; y++) {
    printf("%d ", e_recvbuf[y]);
    currentfield[calcIndex(width, width - 1, y)] = e_recvbuf[y];
  }
  printf("\n");

  printf("Exchanging east bounds %d\n", rank_cart);
  printf("fill east %d ", e_sendbuf);
  for (int y = 0; y < height; y++) {
    e_sendbuf[y] = currentfield[calcIndex(width, width - 2, y)];
    printf("%d ", e_sendbuf[y]);
  }
  printf("\nsend east %d to %d\n", e_sendbuf, w_recvbuf);
  MPI_Sendrecv(e_sendbuf, width, MPI_CHAR, e_id, 2, w_recvbuf, width, MPI_CHAR,
               w_id, 2, cart_comm, MPI_STATUS_IGNORE);
  printf("recv west %d ", w_recvbuf);
  for (int y = 0; y < height; y++) {
    printf("%d ", w_recvbuf[y]);
    currentfield[calcIndex(width, 0, y)] = w_recvbuf[y];
  }
  printf("\n");
}

void game(int width, int height, int num_timesteps, int gsizes[2]) {
  printf("width: %d , height: %d", width, height);
  char *currentfield = calloc(width * height, sizeof(char));
  char *newfield = calloc(width * height, sizeof(char));

  // TODO 1: use your favorite filling
  filling_runner(currentfield, width, height);
  int time = 0;
  write_field(currentfield, gsizes[X], gsizes[Y], time);

  for (time = 1; time <= num_timesteps; time++) {
    printf("\n");
    // TODO 2: implement evolve function (see above)
    evolve(currentfield, newfield, width, height);
    // TODO 3: implement SWAP of the fields
    char *temp = currentfield;
    currentfield = newfield;
    newfield = temp;

    exchange(currentfield, width, height);
    write_field(newfield, gsizes[X], gsizes[Y], time);
  }

  free(currentfield);
  free(newfield);
}

int main(int c, char **v) {
  MPI_Init(&c, &v);

  int width, height, num_timesteps;
  int process_numX;
  int process_numY;

  if (c == 6) {
    width = atoi(v[1]);   ///< read width + 2 boundary cells (low x, high x)
    height = atoi(v[2]);  ///< read height + 2 boundary cells (low y, high y)
    num_timesteps = atoi(v[3]);  ///< read timesteps

    if (width <= 0) {
      width = 32;  ///< default width
    }
    if (height <= 0) {
      height = 32;  ///< default height
    }
    process_numX = atoi(v[4]);  ///< read number of processes in X
    process_numY = atoi(v[5]);  ///< read number of processes in Y

  } else {
    myexit("Too less arguments");
  }

  // TODO get the global rank of the process and save it to rank_global
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);  //-
  // TODO get the number of processes and save it to num_tasks variable
  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);  //-
  printf("MPI processes %d", num_tasks);
  /* Abort if the number of processes does not match with the given
   * configuration.
   */
  if (num_tasks != (process_numX * process_numY)) {
    myexit("ERROR: %d MPI processes needed.\n", process_numX * process_numY);
  }

  /* TODO Create a new cartesian communicator of the worker communicator and get
   * the information.
   */

  int gsizes[2];  // global size of the domain without boundaries
  gsizes[X] = width;
  gsizes[Y] = height;
  int lsizes[2];
  lsizes[X] = width / process_numX;
  lsizes[Y] = height / process_numY;

  int dimsX[2];
  dimsX[X] = process_numX;
  dimsX[Y] = process_numY;
  int periods[2];
  periods[X] = 1;
  periods[Y] = 1;

  MPI_Cart_create(MPI_COMM_WORLD, 2, dimsX, periods, 1, &cart_comm);
  MPI_Comm_rank(cart_comm, &rank_cart);
  MPI_Cart_coords(cart_comm, rank_cart, 2, coords);

  int start_indices_ghost[2];
  start_indices_ghost[X] = coords[X] * lsizes[X];
  start_indices_ghost[Y] = coords[Y] * lsizes[Y];

  /* TODO create and commit a subarray as a new filetype to describe the local
   *      worker field as a part of the global field.
   *      Use the global variable 'filetype'.
   * HINT: use MPI_Type_create_subarray and MPI_Type_commit functions
   */
  int start_indices_file[2];
  start_indices_file[X] = 1;
  start_indices_file[Y] = 1;

  int size_file[2];
  size_file[X] = lsizes[X] - 2;
  size_file[Y] = lsizes[Y] - 2;

  MPI_Type_create_subarray(2, lsizes, size_file, start_indices_file,
                           MPI_ORDER_FORTRAN, MPI_CHAR, &filetype);
  MPI_Type_commit(&filetype);

  /* TODO Create a derived datatype that describes the layout of the inner local
   * field in the memory buffer that includes the ghost layer (local field).
   *      This is another subarray datatype!
   *      Use the global variable 'memtype'.
   */
  MPI_Type_create_subarray(2, gsizes, lsizes, start_indices_ghost,
                           MPI_ORDER_FORTRAN, MPI_CHAR, &memtype);
  MPI_Type_commit(&memtype);

  starts[X] = start_indices_file[X];
  starts[Y] = start_indices_file[Y];

  ends[X] = starts[X] + size_file[X];
  ends[Y] = starts[Y] + size_file[Y];

  printf("starts[x]= %d ,starts[y]= %d \n", starts[X], starts[Y]);
  printf("ends[x]= %d ,ends[y]= %d \n", ends[X], ends[Y]);

  printf("gsizes[x]= %d ,gsizes[y]= %d \n", gsizes[X], gsizes[Y]);
  printf("size_file[x]= %d ,size_file[y]= %d \n", size_file[X], size_file[Y]);

  game(lsizes[X], lsizes[Y], num_timesteps, gsizes);

  MPI_Finalize();
}
