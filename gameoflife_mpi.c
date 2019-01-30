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
// #define CONSOLE_OUTPUT

#define calcIndex(width, x, y) ((y) * (width) + (x))
#define ALIVE 1
#define DEAD 0

#define X 0  //
#define Y 1  //

int rank = 0;       // The current MPI rank in the global communicator.
int rank_cart = 0;  // The current MPI rank in the cart communicator.
int num_tasks;      // The number of processes
int worker_inner_sizes[2];
int worker_sizes[2];
int global_sizes[2];
int coords_of_proc[2];
int process_dim_array[2];
int periodic_boundaries_true[] = {1, 1};  // check

MPI_Datatype filetype;  //
MPI_Comm cart_comm;     // Communicator for the cartesian grid
MPI_File file;          // A shared file pointer
MPI_Datatype memtype;   // A new created type for the inner and outer data
                        // including the ghost layer
MPI_Datatype ghost_layer[4];
MPI_Datatype inner_layer[4];
MPI_Status status;
MPI_Request req;

void myexit(const char* s, ...) {
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
void create_vtk_header(char* header, int width, int height, int timestep) {
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

void write_field(char* currentfield, int width, int height, int timestep) {
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
  MPI_File_delete(filename, MPI_INFO_NULL);

  MPI_File_open(cart_comm, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY,
                MPI_INFO_NULL, &file);
  if (rank_cart == 0) {
    MPI_File_write(file, vtk_header, strlen(vtk_header), MPI_CHAR,
                   MPI_STATUS_IGNORE);
  }
  /* TODO Set the file view for the file handle using collective I/O
   *
   */
  // rc = ...
  int rc;
  rc = MPI_File_set_view(file, header_offset, MPI_CHAR, filetype, "native",
                         MPI_INFO_NULL);
  if (rc != MPI_SUCCESS) {
    myexit("Could not write vtk-Data");
  }
  rc = MPI_File_write_all(file, currentfield, 1, memtype, MPI_STATUS_IGNORE);
  if (rc != MPI_SUCCESS) {
    myexit("write all failed");
  }
  /* TODO Write the data using collective I/O
   *
   */
  // MPI_File_write_all(file, currentfield, sizeof(currentfield) / sizeof(char),
  //                    filetype, &status);
  /* TODO Close the file handle.
   *
   */
  MPI_File_close(&file);
}

void evolve_rank(char* currentfield, char* newfield, int width, int height) {
  // TODO traverse through each voxel and implement game of live logic
  // HINT: avoid boundaries

  int summe_der_Nachbarn;
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      int cell_index = calcIndex(width, x, y);
      newfield[cell_index] = rank;
    }
  }
}
void evolve(char* currentfield, char* newfield, int width, int height) {
  // TODO traverse through each voxel and implement game of live logic
  // HINT: avoid boundaries
  // TODO kommis haben globale starts[]ends[]:
  /*
    starts[X] = start_indices_file[X];
    starts[Y] = start_indices_file[Y];

    ends[X] = starts[X] + size_file[X];
    ends[Y] = starts[Y] + size_file[Y];
    als lokale width und height begrenzung

    trotzdem glsize width und height als index berechnung ??
    int yi = y * width;
    char *nyi = ci + x + (ny * width);
    etc....
    */
  int summe_der_Nachbarn;
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      summe_der_Nachbarn = 0;
      int cell_index = calcIndex(width, x, y);
      // printf("cellindex: %d \n", cell_index);
      // Durchlaufen der 9 Felder des aktuellen "Stempels"
      for (int x1 = -1; x1 <= 1; x1++) {
        for (int y1 = -1; y1 <= 1; y1++) {
          if (currentfield[calcIndex(width, (x + x1), (y + y1))])
            summe_der_Nachbarn++;
          //    printf("summe_der_Nachbarn: %d \n", summe_der_Nachbarn);
        }
      }
      // Wert der untersuchten Zelle von der Summe der Felder abziehen

      if (currentfield[cell_index]) {
        summe_der_Nachbarn--;
      }
      // printf("summe_der_Nachbarn: %d \n", summe_der_Nachbarn);

      // wenn zelle lebt wird 1 von der summe
      // abgezogen

      // if (currentfield[cell_index] == DEAD && summe_der_Nachbarn == 3) {
      //   newfield[cell_index] = ALIVE;
      // }

      if (summe_der_Nachbarn <= 1) {
        newfield[cell_index] = DEAD;
      }

      else if ((summe_der_Nachbarn == 2 && currentfield[cell_index]) == ALIVE ||
               summe_der_Nachbarn == 3) {
        newfield[cell_index] = ALIVE;
      }

      else if (summe_der_Nachbarn >= 4) {
        newfield[cell_index] = DEAD;
      } else {
        newfield[cell_index] = DEAD;
      }
      // HINT: avoid boundaries
    }
  }
}

void filling_random(char* currentfield, int width, int height) {
  int i;
  for (int y = 1; y < height - 1; y++) {
    for (int x = 1; x < width - 1; x++) {
      i = calcIndex(width, x, y);
      currentfield[i] =
          (rand() < RAND_MAX / 10) ? 1 : 0;  ///< init domain randomly
    }
  }
}

void filling_runner(char* currentfield, int width, int height) {
  currentfield[calcIndex(width, width / 2 + 0, height / 2 + 1)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 1, height / 2 + 2)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 2, height / 2 + 0)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 2, height / 2 + 1)] = ALIVE;
  currentfield[calcIndex(width, width / 2 + 2, height / 2 + 2)] = ALIVE;
}
void swap_field(char** currentfield, char** newfield) {
  char* temp = *currentfield;
  *currentfield = *newfield;
  *newfield = temp;
}

void exchange(char* currentfield, int width, int height) {
  int dims[2];
  int periods[2];
  int coords[2];
  int dim;
  int i;
  int r0, r1;  // ranks für sendrcv source und dest.
  int sub_size_up_down[2];
  int sub_size_left_right[2];
  sub_size_up_down[X] = worker_sizes[X];
  sub_size_up_down[Y] = 1;
  sub_size_left_right[X] = 1;
  sub_size_left_right[Y] = worker_sizes[Y];
  int start_ind_ghost[4][2] = {{0, 0},
                               {0, 0},
                               {0, worker_sizes[Y] - 1},
                               {worker_sizes[X] - 1, 0}};  /////////
  int start_ind_inner[4][2] = {{0, 1},
                               {1, 0},
                               {0, worker_sizes[Y] - 2},
                               {worker_sizes[X] - 2, 0}};  // TODO evtl falsch

  // GLayer
  for (int z = 0; z < 4; z++) {
    MPI_Type_create_subarray(1, worker_sizes, sub_size_up_down,
                             start_ind_ghost[z], 1, MPI_CHAR, &ghost_layer[z]);
    MPI_Type_commit(&ghost_layer[z]);
  }
  // inner Layer
  for (int z = 0; z < 4; z++) {
    MPI_Type_create_subarray(1, worker_sizes, sub_size_left_right,
                             start_ind_inner[z], 1, MPI_CHAR, &inner_layer[z]);
    MPI_Type_commit(&inner_layer[z]);
  }
  ///////////
  MPI_Cart_get(cart_comm, 2, dims, periods, coords);
  for (dim = 0, i = 0; dim < dims; ++dim) {
    MPI_Cart_shift(cart_comm, dim, 1, &r0, &r1);
    MPI_Isend(&currentfield, 1, inner_layer[0], r0, MPI_ANY_TAG, cart_comm,
              &req);
    MPI_Irecv(&currentfield, 1, ghost_layer[0], r0, MPI_ANY_TAG, cart_comm,
              &req);

    //  ++i;
    // wieviele req[i] brauche ich?? TODO
    MPI_Isend(&currentfield, 1, inner_layer[2], r1, MPI_ANY_TAG, cart_comm,
              &req);  // eigtl. &req[i]
    MPI_Irecv(&currentfield, 1, ghost_layer[2], r1, MPI_ANY_TAG, cart_comm,
              &req);
    // MPI_Isend(sendbuf + i * sendcount * extent(sendtype), sendcount,
    // sendtype,
    //           r1, ..., cart_comm, &req[i]);
    // MPI_Irecv(recvbuf + i * recvcount * extent(recvtype), recvcount,
    // recvtype,
    //           r1, ..., cart_comm, ...);
    //
    // ++i;

    MPI_Cart_shift(cart_comm, dim, -1, &r0, &r1);
    MPI_Isend(&currentfield, 1, inner_layer[2], r0, MPI_ANY_TAG, cart_comm,
              &req);
    MPI_Irecv(&currentfield, 1, ghost_layer[2], r0, MPI_ANY_TAG, cart_comm,
              &req);

    ++i;  // wieviele req[i] brauche ich?? TODO
    MPI_Isend(&currentfield, 1, inner_layer[0], r1, MPI_ANY_TAG, cart_comm,
              &req);  // eigtl. &req[i]
    MPI_Irecv(&currentfield, 1, ghost_layer[0], r1, MPI_ANY_TAG, cart_comm,
              &req);
    ////////////////////////////////////
    //
    // MPI_Sendrecv(currentfield, 1, inner_layer[0], ghost_layer[2], tag_down,
    //              currentfield, 1, ghost_layer[2], inner_layer[0], tag_up,
    //              cart_comm, MPI_STATUS_IGNORE);
    //              //
    // MPI_Sendrecv(currentfield, 1, inner_layer[2], ghost_layer[0], tag_up,
    //              currentfield, 1, ghost_layer[0], inner_layer[2], tag_down,
    //              cart_comm, MPI_STATUS_IGNORE);
    // // MPI_Sendrecv(array, 1, xr_boundary, nbr_right, tag_right, array, 1,
    // //              xl_boundary, nbr_left, tag_left, MPI_COMM_WORLD,
    // //              MPI_STATUS_IGNORE);
    // // MPI_Sendrecv(array, 1, yd_boundary, nbr_down, tag_down, array, 1,
    // //              yu_boundary, nbr_up, tag_up, MPI_COMM_WORLD,
    // //              MPI_STATUS_IGNORE);

    // /////////////////////////////////////////
  }
  MPI_Waitall(i, req, status);
}
void game(int width, int height, int num_timesteps, int global_sizes[2]) {
  char* currentfield = calloc(width * height, sizeof(char));
  char* newfield = calloc(width * height, sizeof(char));

  // TODO 1: use your favorite filling
  // filling_runner(currentfield, width, height);
  filling_random(currentfield, width, height);
  int time = 0;
  write_field(currentfield, global_sizes[X], global_sizes[Y],
              time);  // TODO kommis hier kommis schreiben globales Feld nicht
                      // worker feld, bshit da game(lsizes) aufruf

  for (time = 1; time <= num_timesteps; time++) {
    // TODO 2: implement evolve function (see above)
    evolve(currentfield, newfield, width,
           height);  // TODO: swap, exchange vor write
    // evolve_rank(currentfield, newfield, width, height);
    swap_field(&currentfield, &newfield);
    exchange(currentfield, width, height);
    write_field(currentfield, global_sizes[X], global_sizes[Y], time);
    // TODO 3: implement SWAP of the fields
  }

  free(currentfield);
  free(newfield);
}

int main(int c, char** v) {
  MPI_Init(&c, &v);
  // if (rank == 0) {
  //   mkdir("./gol/", 0777);
  // }
  int width, height, num_timesteps;  // check
  int process_numX;                  // check
  int process_numY;                  // check
  int worker_size_N_in_x, worker_size_N_in_y;

  int start_worker_inner_offset[2];
  int start_worker_offset[2];
  if (c == 6) {
    process_numX = atoi(v[1]);
    process_numY = atoi(v[2]);
    worker_size_N_in_x = atoi(v[3]);
    worker_size_N_in_y = atoi(v[4]);
    num_timesteps = atoi(v[5]);  ///< read timesteps
    width = worker_size_N_in_x * process_numX;
    // plus zwei für period. RB
    height = worker_size_N_in_y * process_numY;

    if (width <= 0) {
      width = 32;  ///< default width
    }
    if (height <= 0) {
      height = 32;  ///< default height
    }

  } else {
    myexit("Too less arguments");
  }

  // TODO get the global rank of the process and save it  rank_global
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);  //-
  // TODO get the number of processes and save it to num_tasks variable
  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);  //-

  /* Abort if the number of processes does not match with the given
   * configuration.
   */
  if (num_tasks != (process_numX * process_numY)) {
    myexit("ERROR: %d MPI processes needed.\n", process_numX * process_numY);
  }

  /* TODO Create a new cartesian communicator of the worker communicator and
   * get the information.
   */
  process_dim_array[X] = process_numX;
  process_dim_array[Y] = process_numY;

  int ierr = MPI_Cart_create(MPI_COMM_WORLD, 2, process_dim_array,
                             periodic_boundaries_true, 1, &cart_comm);
  // if (MPI_SUCCESS == ierr)
  //   printf("MPI CART Created");
  // else
  //   printf("error");
  MPI_Comm_rank(cart_comm, &rank_cart);
  MPI_Cart_coords(cart_comm, rank_cart, 2, coords_of_proc);

  global_sizes[0] = (worker_size_N_in_x * process_numX);
  global_sizes[1] = (worker_size_N_in_y * process_numY);
  // +2 fuer randaustausch ??notwendig?
  // global size of the domain without boundaries
  worker_inner_sizes[0] = worker_size_N_in_x - 2;
  worker_inner_sizes[1] = worker_size_N_in_y - 2;
  worker_sizes[0] = worker_size_N_in_x;
  worker_sizes[1] = worker_size_N_in_y;
  /* global indices of the first element of the local array */

  start_worker_offset[X] = coords_of_proc[X] * worker_inner_sizes[X];
  // FEHLER bei start... -1
  start_worker_offset[Y] = coords_of_proc[Y] * worker_inner_sizes[Y];

  start_worker_inner_offset[X] = 1;
  // warum absolutes 1 und nicht relatives coords +1
  start_worker_inner_offset[Y] = 1;
  /* TODO create and commit a subarray as a new filetype to describe the local
   *      worker field as a part of the global field.
   *      Use the global variable 'filetype'.
   * HINT: use MPI_Type_create_subarray and MPI_Type_commit functions
   */
  printf(
      "MEMTYPE %d: workersize %d %d, worker_inner_sizes %d %d, start %d %d\n",
      rank, worker_sizes[X], worker_sizes[Y], worker_inner_sizes[X],
      worker_inner_sizes[Y], start_worker_inner_offset[X],
      start_worker_inner_offset[Y]);
  // sleep(100);
  printf(
      "FILETYPE %d: globalsizes %d %d, worker_inner_sizes %d %d, start %d "
      "%d\n",
      rank, global_sizes[X], global_sizes[Y], worker_sizes[X], worker_sizes[Y],
      start_worker_offset[X], start_worker_offset[Y]);
  //========= MEMTYPE =================
  MPI_Type_create_subarray(2, worker_sizes, worker_inner_sizes,
                           start_worker_inner_offset, MPI_ORDER_C, MPI_CHAR,
                           &memtype);  // ERROR
  MPI_Type_commit(&memtype);

  //========  FILETYPE ================
  MPI_Type_create_subarray(2, global_sizes, worker_inner_sizes,
                           start_worker_offset, MPI_ORDER_C, MPI_CHAR,
                           &filetype);
  MPI_Type_commit(&filetype);
  /* TODO Create a derived datatype that describes the layout of the inner
   * local field in the memory buffer that includes the ghost layer (local
   * field). This is another subarray datatype! Use the global variable
   * 'memtype'.
   */

  game(worker_sizes[X], worker_sizes[Y], num_timesteps,
       global_sizes);  // TODO kommis haben worker_sizes[] als w
                       // und h und als 4ten parameter global sizes

  MPI_Finalize();
}
