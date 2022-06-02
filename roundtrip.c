#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define FRACLEVELS 1

int main(int argc, char *argv[])
{
  int mpisize, mpirank;
  int maxFracElems,currentFracLevel;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize); // must be at least 2!! make it work for only one, too? different path?

  // Allocate a 1 MiB buffer
  char *buffer = malloc(sizeof(char) * N);

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
	//buffer = 0,0;
	#pragma omp parallel for num_threads(mpthreads)
	for(int i=maxFracElems; i < mpisize; i++)
	{
		MPI_Send(buffer, 3, MPI_INT, i, i, MPI_COMM_WORLD);
	}
	//print fracElems;
  } else {
       	MPI_Recv(buffer, 3, MPI_INT, 0, mpirank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	//if fraclevel = 0, just ending
	while(currentFracLevel <= FRACLEVELS)
	{
		//calculate number of threads based on mpirank, mpisize, and fraclevel
		//create threads
			//check if new to receive or reuse from previous calculation in this while
			//calculate/proccess
			//send result
		//update currentFracLevel and continue
	}
  }

  MPI_Finalize();
  return 0;
}
