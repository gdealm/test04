#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define FRACLEVELS 2

int main(int argc, char *argv[])
{
  int mpisize, mpirank;
  int maxFracElems,currentFracLevel;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize); // must be at least 2!!

  // Allocate a 1 MiB buffer
  char *buffer = malloc(sizeof(char) * N);

#pragma omp parallel for num_threads(NT) // GGG
for(int i=0; i<NT; i++) printf("thrd no %d of rank %d\n",thrd_no(),mpirank); // GGG

  if (mpirank == 0) {
	int lastLevelElems = pow(2,FRACLEVELS+1);
	int fracElems[lastLevelElems + lastLevelElems - 2]; //(x,y) coordinates for all frac positions
	
	currentFracLevel = 1;
	while(currentFracLevel <= FRACLEVELS)
	{
		int mpthreads = currentFracLevel*2;
		#pragma omp parallel for num_threads(mpthreads)
		for(int i=0; i < mpthreads; i++)
		{
			MPI_INT[] buffer = {currentFracLevel,xPos[],yPos[]};
			int mprank = omp_get_thread_num();
			if(mprank % 2 = 0)
			{
				if(mprank <= maxFracElems)
				{
					MPI_Send(buffer, 3, MPI_INT, ((((mprank/2)+1)%(mpisize-1))+1), ((mprank/2)+1), MPI_COMM_WORLD);
				}
				MPI_Recv(buffer, 2, MPI_INT, ((((mprank/2)+1)%(mpisize-1))+1), (((mprank/2)+1)*10), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			else
			{
				MPI_Recv(buffer, 2, MPI_INT, ((((mprank/2)+1)%(mpisize-1))+1), ((((mprank/2)+1)*10)+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
		prevFracLevel = currentFracLevel;
		currentFracLevel++;
	}
	//activate not activated yet, if any, to start to finish 

  } else {
       	MPI_Recv(buffer, 3, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	omp_get_thread_num();
	/*
        printf("Rank %d (running on '%s'): receive the message and sending it to rank %d\n",rank,hostname,(rank+1)%size);
	MPI_Send(buffer, N, MPI_BYTE, (rank+1)%size, 1, MPI_COMM_WORLD);
	*/
  }

  MPI_Finalize();
  return 0;
}
