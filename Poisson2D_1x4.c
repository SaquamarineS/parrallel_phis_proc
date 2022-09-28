#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#define N 50
// #define Ny 100
#include <math.h>
#include <string.h>


int main(int argc, char *argv[]) {
    double t1 = MPI_Wtime();
    int Ny = N*2;
    int rank, size;
    double dx = 1.0/(N*4);
    double dy = 1.0/(Ny);
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;

    // store some values in x and y here

    double(*locgrid)[N] = malloc (sizeof(double[Ny][N]));

    if(locgrid == NULL)
    {
        printf("ERRRORR MALLOC\n");
    }
    for (int i =1; i<Ny-1; i++){
        for (int j = 1; j<N-1; j++){
            locgrid[i][j] = 0; // IV fill
        }
    }

    double(*nextlocgrid)[N] = malloc (sizeof(double[Ny][N]));
    if(nextlocgrid == NULL)
    {
        printf("ERRRORR MALLOC\n");
    }
    for (int i =1; i<Ny-1; i++){
        for (int j = 1; j<N-1; j++){
            nextlocgrid[i][j] = 0; // IV fill
        }
    }
    double varepsilon = 0.001;
    int k = 0; 
    while(1){
        k++;
        // boundaries

        // 
        if (rank ==0){
            for (int i = 0; i<=.3*Ny; i++){ // heat 1 left
                locgrid[i][0] = 1;
            }
            for (int i = 0.6*Ny; i<=.8*Ny; i++){ // heat 2left
                locgrid[i][0] = 1;
            }
            // for (int i = 0.6*N; i<N; i++){ // heat bot
            //     locgrid[0][i] = 1;
            // } 
            
        }
        if (rank == 1){
            for (int i = 0.2*N; i<N; i++){ // heat bot
                locgrid[0][i] = 1;
            } 
            // printf("rank%d: %d\n",rank,k);
            // for (int i =0; i<Ny; i++){
            //     for (int j = 0; j<N; j++){
            //         printf("%4.1f ", locgrid[i][j]);
            //     }
            //     printf("\n");
            // }  
        }

        if (rank == 2){
            for (int i = 0; i<0.4*N; i++){ // heat bot
                locgrid[0][i] = 1;
            }
            // printf("rank%d: %d\n",rank,k);
            // for (int i =0; i<Ny; i++){
            //     for (int j = 0; j<N; j++){
            //         printf("%4.1f ", locgrid[i][j]);
            //     }
            //     printf("\n");
            // }
        }
        if (rank == 3){
                    // for (int i = 0; i<0.2*N; i++){
                    //     locgrid[0][i] = 1;                  // heater bot
                    // }
            for (int i = 0; i<0.2*Ny; i++){
                locgrid[i][N-1] = locgrid[i][N-2]; // hole 1 right
            }
            for (int i = 0.4*Ny; i<0.6*Ny; i++){
                locgrid[i][N-1] = locgrid[i][N-2]; // hole 2 right
            }
            for (int i = 0.7*Ny; i<0.9*Ny; i++){
                locgrid[i][N-1] = locgrid[i][N-2]; // hole 3 right
            }
            // printf("rank%d: %d\n",rank,k);
            // for (int i =0; i<Ny; i++){
            //     for (int j = 0; j<N; j++){
            //         printf("%4.1f ", locgrid[i][j]);
            //     }
            //     printf("\n");
            // }
        }
        for (int i = 0; i < Ny; i++){
            for (int j = 0; j < N; j++){
                nextlocgrid[i][j]=locgrid[i][j];
            }
        }
        // solving local poisson
        int xstart = N * rank;
        for (int i = 1; i < Ny-1; i++){
            for (int j = 1; j < N-1; j++){
                // nextlocgrid[i][j] = 0.25 * (locgrid[i+1][j]+locgrid[i-1][j]+locgrid[i][j+1]+locgrid[i][j-1]);
                nextlocgrid[i][j] = ((locgrid[i+1][j]+locgrid[i-1][j])/dy/dy+(locgrid[i][j+1]+locgrid[i][j-1])/dx/dx+((xstart + j)*dx)*((xstart + j)*dx))/(2/dx/dx+2/dy/dy);
            }
        }

        // printf("time for communication\n");
        if (rank == 0){
            double *commbuffer;
            commbuffer = malloc(Ny*sizeof(double));
            for (int i = 0; i<Ny-2; i++){
                commbuffer[i] = locgrid[i+1][N-2];
            }
            // printf("my rank: %d, i m sending  to neighb", rank);
            MPI_Rsend(commbuffer, Ny-1, MPI_DOUBLE,1,1234,MPI_COMM_WORLD);
            
            MPI_Recv(commbuffer,Ny-1,MPI_DOUBLE, 1,1234,MPI_COMM_WORLD,&status);
            for (int i = 0; i<Ny-2; i++){
                locgrid[i+1][N-1] = commbuffer[i];
            }   
        }
        if (rank == 1){
            double *commbuffer;
            commbuffer = malloc(Ny*sizeof(double));
            MPI_Recv(commbuffer,Ny-1,MPI_DOUBLE, 0,1234,MPI_COMM_WORLD,&status);
            for (int i = 0; i<Ny-2; i++){
                locgrid[i+1][0] = commbuffer[i];
            }
            for (int i = 0; i<Ny-2; i++){
                commbuffer[i] = locgrid[i+1][1];
            }
            MPI_Rsend(commbuffer, Ny-1, MPI_DOUBLE,0,1234,MPI_COMM_WORLD);




            MPI_Recv(commbuffer,Ny-1,MPI_DOUBLE, 2,1234,MPI_COMM_WORLD,&status);
            for (int i = 0; i<Ny-2; i++){
                locgrid[i+1][N-1] = commbuffer[i];
            }
            for (int i = 0; i<Ny-2; i++){
                commbuffer[i] = locgrid[i+1][N-2];
            }
            MPI_Rsend(commbuffer, Ny-1, MPI_DOUBLE,2,1234,MPI_COMM_WORLD);

        }
        if (rank == 2){
            double *commbuffer;
            commbuffer = malloc(Ny*sizeof(double));

            for (int i = 0; i<Ny-2; i++){
                commbuffer[i] = locgrid[i+1][N-2];
            }
            MPI_Rsend(commbuffer, Ny-1, MPI_DOUBLE,3,1234,MPI_COMM_WORLD);
            MPI_Recv(commbuffer,Ny-1,MPI_DOUBLE, 3,1234,MPI_COMM_WORLD,&status);
            for (int i = 0; i<Ny-2; i++){
                locgrid[i+1][N-1] = commbuffer[i];
            }

            for (int i = 0; i<Ny-2; i++){
                commbuffer[i] = locgrid[i+1][1];
            }
            MPI_Rsend(commbuffer, Ny-1, MPI_DOUBLE,1,1234,MPI_COMM_WORLD);

            MPI_Recv(commbuffer,Ny-1,MPI_DOUBLE, 1,1234,MPI_COMM_WORLD,&status);
            for (int i = 0; i<Ny-2; i++){
                locgrid[i+1][0] = commbuffer[i];
            }

            
            
        }

        if (rank == 3){
            double *commbuffer;
            commbuffer = malloc(Ny*sizeof(double));
            MPI_Recv(commbuffer,Ny-1,MPI_DOUBLE, 2,1234,MPI_COMM_WORLD,&status);
            for (int i = 0; i<Ny-2; i++){
                locgrid[i+1][0] = commbuffer[i];
            }
            for (int i = 0; i<Ny-2; i++){
                commbuffer[i] = locgrid[i+1][1];
            }
            MPI_Rsend(commbuffer, Ny-1, MPI_DOUBLE,2,1234,MPI_COMM_WORLD);
        }


        // double *maxDelta;
        // maxDelta = malloc(size * sizeof(double));
        double maxDelta = fabs(locgrid[0][0] - nextlocgrid[0][0]);
        // if (rank == 1){
        //     printf("maxdelta changed to: %f\n ",nextlocgrid[0][0]);
        // }
        for (int i = 0; i < Ny; i++){
            for (int j = 0; j < N; j++){
                if (fabs(locgrid[i][j] - nextlocgrid[i][j])>maxDelta){
                    maxDelta = fabs(locgrid[i][j] - nextlocgrid[i][j]);
                    // printf("maxdelta changed to: %f\n ",maxDelta);
                }
                // printf("maxdelta changed to: %f\n ",fabs(locgrid[i][j] - nextlocgrid[i][j]));
            }
        }
        // gathering local maxes
        if (rank == 1){
            MPI_Rsend(&maxDelta, 1, MPI_DOUBLE,0,1234,MPI_COMM_WORLD);
            // printf("%d: %d iter, maxdelta = %f\n", rank, k, maxDelta);
        }

        if (rank == 2){
            MPI_Rsend(&maxDelta, 1, MPI_DOUBLE,0,1234,MPI_COMM_WORLD);
            // printf("%d: %d iter, maxdelta = %f\n", rank, k, maxDelta);
        }
        if (rank == 3){
            MPI_Rsend(&maxDelta, 1, MPI_DOUBLE,0,1234,MPI_COMM_WORLD);
            // printf("%d: %d iter, maxdelta = %f\n", rank, k, maxDelta);
        }
        double maxDelta2, maxDelta3, maxDelta4;
        int breakflag=0;
        if (rank ==0 ){
            for (int i = 1; i<size; i++){
                MPI_Recv(&maxDelta2,1,MPI_DOUBLE, i,1234,MPI_COMM_WORLD,&status);
                if (maxDelta>maxDelta2){
                    if (maxDelta < varepsilon){
                        breakflag = 1;
                        // printf("%f \n", maxDelta);
                    }
                }
                if (maxDelta2 > maxDelta){
                    maxDelta = maxDelta2;
                    if (maxDelta < varepsilon){
                        breakflag = 1;
                        // printf("%f \n", maxDelta);
                    }
                }
            }
            // if (k%1000==0){
            //     // printf("%d iter, maxdelta = %f\n", k, maxDelta);
            // }
        }
        MPI_Bcast(&breakflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
        for (int i = 1; i < Ny-1; i++){
            for (int j = 1; j < N-1; j++){
                locgrid[i][j] = nextlocgrid[i][j];
            }
        }
        // printf("%d , rank: %d\n",breakflag, rank);
        if (breakflag==1){
            if (rank == 0){
                double t2 = MPI_Wtime();
                printf("my rank is: %d\nPoisson eq'n CONVERGED with Jacobi method\nunder varepsilon = %f, maxerror = %f \ntime taken = %12.8f seconds\n", 
                rank,varepsilon, maxDelta, t2-t1);
                printf("writing output files to \"/Users/user/Documents/MATLAB/\"\n");
            }
            char* a = "/Users/user/Documents/MATLAB/OUTnp4_dim1_"; 
            char* extension =".csv";
            char fileSpec[strlen(a)+strlen("_rank")+1+strlen(extension)+1];
            FILE *out;
            snprintf( fileSpec, sizeof( fileSpec ), "%s_rank%d%s", a,rank, extension );
            out = fopen( fileSpec, "w+" );


            for(int i = 1; i<Ny-2; i ++){
                for (int j = 1; j<N-2; j++){
                    fprintf(out,"%7.5f,",locgrid[i][j]);
                }
                fprintf(out,"\n");
            }
            fclose(out);
            
            free(locgrid);
            free(nextlocgrid);
            MPI_Finalize();
            return 0;
        }

    }

    free(locgrid);
    MPI_Finalize(); 
    return 0;
}
