#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
// 2x2 topology, cyclic in x direction
# define INNER_N 2 // 0 to 9 total of 10


int main(int argc, char *argv[]){
    int rank, size;
    MPI_Comm comm1;
    int outer_n = INNER_N + 2;
    int dimvec[2], periodvec[2], reorder;
    int id;
    int dimensions;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Status Status;
    // if (size !=6){
    //     printf("please run with 4 processess.\n"); fflush(stdout);
    //     MPI_Abort(MPI_COMM_WORLD,1);
    // }
    dimensions = 2;
    dimvec[0] = 2; dimvec[1] = 2;
    periodvec[0] = 0; periodvec[1] = 0; // tut nnado 0 postavit! net perioda
    reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD,dimensions,dimvec,periodvec, reorder, &comm1);
    

    // if (rank==0){
    //         printf("\n 2x2 comm1\n");
    // }
    int coord[2];
    int findrank;
    int findcoord;
    int maxdims = 2;
    for (findcoord =0; findcoord <size; findcoord++){
        MPI_Cart_coords(comm1, findcoord, maxdims, coord);
        MPI_Cart_rank(comm1, coord, &findrank);
        // if (rank == 0){
        //     printf("for rank:%d found coord:(%d; %d); for coord found rank:%d\n", findcoord, coord[0],coord[1], findrank);
        // } 
        if (rank==findcoord){
            int up, down, right, left;
            MPI_Cart_shift(comm1, 0, 1, &left, &right);
            MPI_Cart_shift(comm1, 1, 1, &up, &down);
            // printf("%d: up=%d; down=%d; left=%d; right=%d\n", rank, up, down, left, right);
        }
    }
    
    // if (rank==0){
    //         printf("\n 1x4 comm2\n");
    // }
    MPI_Comm comm2;
    dimensions = 2;
    dimvec[0] = 4; dimvec[1] = 1;
    periodvec[0] = 0; periodvec[1] = 0; // tut nnado 0 postavit! net perioda
    reorder = 1;
    MPI_Cart_create(MPI_COMM_WORLD,dimensions,dimvec,periodvec, reorder, &comm2);
    maxdims = 2;

    for (findcoord =0; findcoord <size; findcoord++){
        MPI_Cart_coords(comm2, findcoord, maxdims, coord);
        MPI_Cart_rank(comm2, coord, &findrank);
        // if (rank == 0){
        //     printf("for rank:%d found coord:(%d; %d); for coord found rank:%d\n", findcoord, coord[0],coord[1], findrank);
        // } 
        // if (rank==findcoord){
        //     int up, down, right, left;
        //     MPI_Cart_shift(comm1, 0, 1, &left, &right);
        //     MPI_Cart_shift(comm1, 1, 1, &up, &down);
        //     printf("%d: up=%d; down=%d; left=%d; right=%d\n", rank, up, down, left, right);
        // }
    }
    
    // boundaries: 

    // double *localgrid;
    // localgrid = malloc(size * INNER_N * INNER_N * sizeof(double));

    double localgrid[outer_n][outer_n];
    int tag = 1;
    int tag2 = 2;
    if (size==4){
        // IV
        // fill boundary 


        if (rank == 0){
            for (int i=1; i < outer_n; i++){
                for (int j = 1; j < outer_n; j++){
                    localgrid[i][j] = i+j*2; // fill with IVal 
                }
            }
            // printout local grid before comm
            for (int i=0; i <= outer_n; i++){
                for (int j = 0; j <= outer_n; j++){
                    printf("%4.1f ",localgrid[i][j]);
                }
                printf("\n");
            }

            // SOLVE LOCAL POISSON HERE


            //ghost right - to proc 2
            double sendright[INNER_N];
            for (int i = 1; i< outer_n; i++){
                sendright[i-1] = localgrid[i][outer_n - 1];
                // printf("[%2d] = %4.1f ",i,sendright[i-1]);
            }
            printf("\n");
            //ghost down - to proc 1
            double senddown[INNER_N];
            for (int i = 1; i< outer_n; i++){
                senddown[i-1] = localgrid[outer_n - 1][i];
                // printf("[%2d] = %4.1f ",i,senddown[i-1]);
            }
            
            int neighX = 2; 
            int neighY = 1; 
            double *arrXrecv;
            arrXrecv = malloc(INNER_N * sizeof(double));
            double arrYrecv[INNER_N];
            MPI_Rsend(sendright, INNER_N+1, MPI_DOUBLE,neighX,tag,MPI_COMM_WORLD);
            MPI_Recv(arrXrecv,INNER_N+1,MPI_DOUBLE, neighX,tag2,MPI_COMM_WORLD,&Status);
            for (int i = 1; i < outer_n; i++){
                localgrid[i][outer_n] = arrXrecv[i-1];
            }
            for (int i = 0; i < 3; i++){
                printf("%f \n", arrXrecv[i]);
            }
            // y communication
            // for (int i = 1; i < outer_n; i++){
            //     localgrid[outer_n][i] = arrYrecv[i-1];
            // }

            printf("rank0\n");
            for (int i=0; i <= outer_n; i++){
                for (int j = 0; j <= outer_n; j++){
                    printf("%4.1f ",localgrid[i][j]);
                }
                printf("\n");
            }
            
        }
        if (rank == 2){
            printf("rank2\n");
            for (int i=1; i < outer_n; i++){
                for (int j = 1; j < outer_n; j++){
                    localgrid[i][j] = i*2+j; // fill with IVal 
                }
            }
            // printout local grid
            for (int i=0; i <= outer_n; i++){
                for (int j = 0; j <= outer_n; j++){
                    printf("%4.1f ",localgrid[i][j]);
                }
                printf("\n");
            }

            // SOLVE LOCAL POISSON HERE

            //ghost right - to proc 2
            double *sendleft;
            sendleft = malloc(INNER_N * sizeof(double));
            for (int i = 0; i<= INNER_N; i++){
                sendleft[i] = localgrid[i+1][1];
                // printf("%4.1f ",localgrid[i+1][1]);
            }
            printf("\n");
            //ghost down - to proc 1
            double senddown[INNER_N];
            for (int i = 1; i< outer_n; i++){
                senddown[i-1] = localgrid[outer_n - 1][i];
                // printf("[%2d] = %4.1f ",i,senddown[i-1]);
            }
            
            double arrXrecv[INNER_N];
            double arrYrecv[INNER_N];
            int neighX = 0; 
            int neighY = 3;
            // int to

            MPI_Recv(arrXrecv,INNER_N+1,MPI_DOUBLE, neighX,tag,MPI_COMM_WORLD,&Status);
            MPI_Rsend(sendleft, INNER_N+1, MPI_DOUBLE,neighX,tag2,MPI_COMM_WORLD);
            free(sendleft);
        }

        if (rank == 1){
            
        }
        
        if (rank == 3){
            
        }
    }


    MPI_Finalize();

} 
