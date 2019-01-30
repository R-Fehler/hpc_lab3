#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>

#define calcIndex(width, x,y)  ((y)*(width) + (x))

#define EPS 0.25
#define EPS_R (1.0/0.25)
#define TAU 1.0
#define TAU_R (1.0/TAU)
#define dt 0.0001
#define dx 0.1
#define dx_R2 (1.0/(dx*dx))
#define gamma 1.0
#define DRIVING_FORCE -90

//MPI 
int rank;
int procs;
int ndims = 1;
int dims[1];
int periods[1];
int rank_left;
int rank_right;  
int coords[1];
MPI_Comm   comm_cart;
MPI_Status status; 


MPI_Datatype ghostlayerleft; 
MPI_Datatype ghostlayerright; 
MPI_Datatype innerlayerleft; 
MPI_Datatype innerlayerright; 


void writeVTK2(long timestep, double* data, char prefix[1024], long w, long h) {
  char filename[2048];  
  long x, y; 
  
  long offsetX = coords[0]*(w-3);
  long offsetY = 0;
  float deltax = 1.0;
  float deltay = 1.0;
  long  nxy = w * h * sizeof(float);  

  snprintf(filename, sizeof(filename), "%s%d-%05ld%s", prefix, coords[0], timestep, ".vti");
  FILE* fp = fopen(filename, "w");

  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp, "<ImageData WholeExtent=\"%ld %ld %ld %ld %ld %ld\" Origin=\"0 0 0\" Spacing=\"%le %le %le\">\n", offsetX, offsetX + w, offsetY, offsetY + h, (long)0, (long)0, deltax, deltay, 0.0);
  fprintf(fp, "<CellData Scalars=\"%s\">\n", prefix);
  fprintf(fp, "<DataArray type=\"Float32\" Name=\"%s\" format=\"appended\" offset=\"0\"/>\n", prefix);
  fprintf(fp, "</CellData>\n");
  fprintf(fp, "</ImageData>\n"); 
  fprintf(fp, "<AppendedData encoding=\"raw\">\n");
  fprintf(fp, "_");
  fwrite((unsigned char*)&nxy, sizeof(long), 1, fp);

  for (y = 0; y < h; y++) {
    for (x = 0; x < w; x++) {
      float value = data[calcIndex(w, x,y)];
      fwrite((unsigned char*)&value, sizeof(float), 1, fp);
    }
  }
  
  fprintf(fp, "\n</AppendedData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}


void show(double* currentfield, long w, long h) {
  printf("\033[H");
  int x,y;
  for (y = 0; y < h; y++) {
    for (x = 0; x < w; x++) printf(currentfield[calcIndex(w, x,y)] > 0.5 ? "\033[07m  \033[m" : "  ");
    //printf("\033[E");
    printf("\n");
  }
fflush(stdout);
  usleep(SLEEPTIME);
}

void filling(double* currentfield, long w, long h) {
  for (long y = 0 /*h/2-h/8*/; y < /*h/2 +*/  h/4; y++) {
    for (long x = 0 /*w/2-w/8*/; x < /*w/2 + */ w/4; x++) {
      currentfield[calcIndex(w, x,y)] = 1.0;
    }
  }
}
 
void phasefield_sweep(double* field_src, double* field_dst, long w, long h) {
  for (long y = 1; y < h - 1; y++) {
    for (long x = 1; x < w - 1; x++) {
//       IACA_START
      long index = calcIndex(w, x,y);
      field_dst[index] = field_src[index] + dt*TAU_R * 
      (
        gamma * 
        (  
            (field_src[index+1] - 2*field_src[index] +  field_src[index-1])*(dx_R2) +  
            (field_src[index+w] - 2*field_src[index] +  field_src[index-w])*(dx_R2) 
        )
        -  gamma * (2.0 * field_src[index] * field_src[index] * field_src[index] -3.0 * field_src[index] * field_src[index] + field_src[index]) 
        -  DRIVING_FORCE * field_src[index] * (1.0 - field_src[index])
      );
//       IACA_END
    }
  }
}

void copycells(double* field, long w, long h, long fromx, long fromy, long tox, long toy, long width, long height) {
  for (long y = 0; y < height; y++) {
    for (long x = 0; x < width; x++) {
      field[calcIndex(w, tox+x,toy+y)] = field[calcIndex(w, fromx+x,fromy+y)];
    }
  }
}

void ghostlayer_exchange(double* field) {
  MPI_Request request[4];
  MPI_Status  status[4];
  MPI_Irecv(field, 1, ghostlayerright, rank_right, 1, comm_cart, &(request[0]));     
  MPI_Irecv(field, 1, ghostlayerleft,  rank_left,  2, comm_cart, &(request[1]));
  MPI_Isend(field, 1, innerlayerleft,  rank_left,  1, comm_cart, &(request[2]));
  MPI_Isend(field, 1, innerlayerright, rank_right, 2, comm_cart, &(request[3]));
  
  MPI_Waitall(4, request, status); 
}

// void ghostlayer_exchange(double* field) {
//   if (coords[0] % 2 == 0) {
//         MPI_Send(field, 1, innerlayerright, rank_right, 1, comm_cart);
//         MPI_Recv(field, 1, ghostlayerright, rank_right, 2, comm_cart, &status);     
//         MPI_Send(field, 1, innerlayerleft, rank_left, 3, comm_cart);
//         MPI_Recv(field, 1, ghostlayerleft, rank_left, 4, comm_cart, &status);
//     } else {
//         MPI_Recv(field, 1, ghostlayerleft, rank_left, 1, comm_cart, &status);
//         MPI_Send(field, 1, innerlayerleft, rank_left, 2, comm_cart);
//         MPI_Recv(field, 1, ghostlayerright, rank_right, 3, comm_cart, &status);
//         MPI_Send(field, 1, innerlayerright, rank_right, 4, comm_cart);
//     }
// }

void boundarycondition_periodic(double* field, long w, long h) {
   //            domain   from      to      size
   copycells(field, w,h,  0,  h-2,  0,  0,  w,1); // top to buttom
   copycells(field, w,h,  0,  1,    0,h-1,  w,1); // buttom to top
}

void initMPI(long w, long h) {
   //init
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Identify this process
  MPI_Comm_size(MPI_COMM_WORLD, &procs); // Find out how many total processes are active
  
  periods[0] = 1;     //numner of dimensions
  dims[0]    = procs; //number of processes in x dim

  MPI_Cart_create(
    MPI_COMM_WORLD, //[in] input communicator
    ndims,          //[in] number of dimensions of cartesian grid
    dims,           //[in] integer array of size ndims specifying the number of processes in each dimension 
    periods,        //[in] logical array of size ndims specifying whether the grid is periodic (true) or not (false) in each dimension 
    1,              //[in] ranking may be reordered (true) or not (false)
    &comm_cart      //[out] communicator with new cartesian topology 
  );
  
  // get neighbor ranks
  MPI_Cart_shift(
    comm_cart,  //[in] communicator with cartesian structure
    0,          //[in] coordinate dimension of shift
    1,          //[in] displacement (> 0: upwards shift, < 0: downwards shift)
    &rank_left, //[out] rank of source process
    &rank_right  //[out] rank of destination process 
  );
    
  MPI_Cart_coords(
    comm_cart, //[in] communicator with cartesian structure
    rank,      //[in] rank of a process within group of comm
    ndims,   //[in] length of vector coords in the calling program
    coords    //[out] integer array (of size ndims) containing the Cartesian coordinates of specified process
  );
  
//   printf("rank %d, coord %d\n", rank, coords[0]);
  int gsizes[2] = {(int)w, (int)h};
  int lsizes[2] = {1,      (int)h};
  
  int start_indices[2] = {0, 0};
  MPI_Type_create_subarray(2, gsizes, lsizes, start_indices, MPI_ORDER_FORTRAN, MPI_DOUBLE, &ghostlayerleft);
  MPI_Type_commit(&ghostlayerleft); 
  
  start_indices[0] = 1;
  MPI_Type_create_subarray(2, gsizes, lsizes, start_indices, MPI_ORDER_FORTRAN, MPI_DOUBLE, &innerlayerleft);
  MPI_Type_commit(&innerlayerleft); 
  
  start_indices[0] = w-1;
  MPI_Type_create_subarray(2, gsizes, lsizes, start_indices, MPI_ORDER_FORTRAN, MPI_DOUBLE, &ghostlayerright);
  MPI_Type_commit(&ghostlayerright); 
  
  start_indices[0] = w-2;
  MPI_Type_create_subarray(2, gsizes, lsizes, start_indices, MPI_ORDER_FORTRAN, MPI_DOUBLE, &innerlayerright);
  MPI_Type_commit(&innerlayerright); 
  
}


void solverloop(long w, long h, long timesteps) {

  double *field_src = calloc(w*h, sizeof(double));
  double *field_dst = calloc(w*h, sizeof(double));
  
  filling(field_src, w, h);
  boundarycondition_periodic(field_src, w, h);
  ghostlayer_exchange(field_src);
  writeVTK2(-1,field_src,"phi1", w, h);

  for (long t=0; t<timesteps; t++) {
    phasefield_sweep(field_src, field_dst, w, h);
    boundarycondition_periodic(field_dst, w, h);
    ghostlayer_exchange(field_dst);
    
//    SHOWFIELD(show(field_dst, w, h);)
    //TODO secure by tokening for >>N
    WRITTE_VTK(writeVTK2(t,field_src,"heat", w, h);)
    PRINT_TIMEITER(printf("time %f at %ld timestep\n",t*dt, t);)

    double *temp = field_src;     //SWAP
    field_src = field_dst;
    field_dst = temp;
  }
  //deinit
  free(field_src);
  free(field_dst);
}
 
int main(int c, char **v) {
  /* Initialize the infrastructure necessary for communication */
  MPI_Init(&c, &v);
  long w = 20, h = 20, t = 100; //defaults
  if (c > 1) w = atoi(v[1]); ///< read width
  if (c > 2) h = atoi(v[2]); ///< read height
  if (c > 3) t = atoi(v[3]); ///< read time steps
  initMPI(w, h);
  solverloop(w, h, t);
  
  MPI_Finalize();
}
