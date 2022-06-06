#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <mpi.h>
#include <omp.h>

#define FRACLEVELS 1 // Define here the number of levels the fractal tree will have

int mypow(int base, int exp)
{
	int result = 1;
	for(int i = 1; i <= exp; i++)
	{
		result *= base;
	}
	return result;
}

int main(int argc, char *argv[])
{
  int mpisize; // number of MPI processing machines
  int mpirank; // MPI rank of this machine
	
  int currFracLevel; // current fractal level being treated
  int maxFracElems = 0; // max number of fractal elements processed so far in parallel

  MPI_Init(&argc, &argv); // initialize MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank); // obtain mpirank
  MPI_Comm_size(MPI_COMM_WORLD, &mpisize); // obtain mpisize //GGG if mpisize = 1, it will not work

  printf("Rank %d of %d is starting. \n",mpirank,mpisize);	
	
  // controller task

  if (mpirank == 0) {
	int lastLevelElems = mypow(2,FRACLEVELS-1); // calculate the number of elements in the last level of the fractal
	int fracElems[lastLevelElems + lastLevelElems - 1][2]; // define array to contain (x,y) coordinates of all the fractal elements
	int padval = 5; // value to calculate padding from the borders for the elements
	
	currFracLevel = 1; // just to indicate it will be calculating the level 1 single element positions
	fracElems[0][0] = 0; // x position
	fracElems[0][1] = lastLevelElems-1; // y position
	
	currFracLevel = 2; // indicates current processing level 2 to start the loop
	/*
	// starts calculating levels from level 2 to the last one defined
	while(currFracLevel <= FRACLEVELS)
	{
		int mpthreads = mypow(2,currFracLevel-1); // number of threads for the current level is the number of elements in the level / 2
		#pragma omp parallel for num_threads(mpthreads)
		for(int i=0; i < mpthreads; i++)
		{
			int buffer[2]; // create a buffer to receive calculated position
			int mprank = omp_get_thread_num(); // OpenMP rank of this thread
			//if mprank is higher than last level processed, a new "virtual" MPI machine will start processing some branch
			if(mprank >= maxFracElems)
			{
				int sendBuffer[3];  // buffer to send: current level, parent element x position and parent element y position
				int posOrigin = ((mpthreads/2)-1) + (mprank/2); // calculate index of the parent element
				sendBuffer[0] = currFracLevel; // set current level 
				sendBuffer[1] = fracElems[posOrigin][0]; // set parent element x position
				sendBuffer[2] = fracElems[posOrigin][1]; // set parent element y position
				//MPI_Send(sendBuffer, 3, MPI_INT, (((mprank/2)%(mpisize-1))+1), ((mprank/2)+1), MPI_COMM_WORLD); 
				MPI_Send(sendBuffer, 3, MPI_INT, ((mprank%(mpisize-1))+1), (mprank+1), MPI_COMM_WORLD); // send to next MPI machine in round robin
			}
			MPI_Recv(buffer, 2, MPI_INT, ((mprank%(mpisize-1))+1), (mprank+1), MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive calculated element position
			fracElems[mpthreads - 1 + mprank][0] = buffer[0]; // update element x position in consolidated array
			fracElems[mpthreads - 1 + mprank][1] = buffer[1]; // update element y position in consolidated array
		}
		maxFracElems = mpthreads; // update max threads for next level control
		currFracLevel++; // increment frac level to the next one
	}
	*/
	//currFracLevel = FRACLEVELS + 1; //GGG remove
	//activate not activated yet, if any, to continue from receive to end
	if(maxFracElems < mpisize-1)
	{
		/*
		int sendBuffer[3];  // buffer to send: current level
		sendBuffer[0] = currFracLevel; // current fractal level is above to be treated already
		sendBuffer[1] = 0; // just to initialize
		sendBuffer[2] = 0; // just to initialize
		*/
		//send buffer to MPI machines not activated yet
		#pragma omp parallel for num_threads((mpisize-1)-maxFracElems)
		for(int i=maxFracElems; i < mpisize-1; i++)
		{	
			printf("send %d\n",(i+1));
			int sendBuffer[3];  // buffer to send: current level
			sendBuffer[0] = currFracLevel; // current fractal level is above to be treated already
			sendBuffer[1] = 0; // just to initialize
			sendBuffer[2] = 0; // just to initialize
			MPI_Send(sendBuffer, 3, MPI_INT, i+1, i+1, MPI_COMM_WORLD);
		}
		printf("all sent\n");	
	}
	// print fractal elements positions;
	int fracElemsLength = sizeof fracElems / sizeof fracElems[0];
	for(int i = 0; i < fracElemsLength; i++)
	{
		printf("(%d,%d)\n",fracElems[i][0],fracElems[i][1]);	
	}
  }  else { // not MPI member 0, calculate element(s) to send to 0
	int buffer[3]; // create a buffer to receive parent element position
	printf("receive %d\n",mpirank);
       	MPI_Recv(buffer, 3, MPI_INT, 0, mpirank, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive first message from MPI member 0
	printf("received %d\n",mpirank);
	currFracLevel = buffer[0]; // check in which level was first activated
	/*
	if(currFracLevel <= FRACLEVELS) // will do something if still in some level to process
	{
		bool first = true; // used to indicate if it is the first level being treated or not to treat the first element as a new receival further
		int displacementX = ((FRACLEVELS-currFracLevel)+1)*2; // to add or subtract from parent x position to calculate the x position for this level elements
		int fracLevelElems = mypow(2,FRACLEVELS-1); // calculate the number of elements in the last level of the fractal
		int lastLevelLocalElems = (fracLevelElems/(mpisize-1))+1; // calculate +- the number of elements in the last level of the fractal this MPI will treat
		int localFracElems[lastLevelLocalElems][2];
		localFracElems[0][0] = buffer[1];
		localFracElems[0][1] = buffer[2];
		maxFracElems = 1;
		while(currFracLevel <= FRACLEVELS)
		{
			fracLevelElems = mypow(2,currFracLevel-1); // calculate the number of elements in the current level of the fractal
			int mpthreads = (fracLevelElems/(mpisize-1)); // contains the number of elements in the current level of the fractal this MPI will treat
			if(mpirank < (fracLevelElems%(mpisize-1))+1) // may need to add one more if this MPI is treatment some of the remainder
			{	
				mpthreads++;
			}
			//create threads
			#pragma omp parallel for num_threads(mpthreads)
			for(int i=0; i < mpthreads; i++)
			{
				int mprank = omp_get_thread_num(); // OpenMP rank of this thread
				//check if new to receive or reuse from previous calculation in this while
				if(mprank >= maxFracElems)
				{
					int buffer[3]; // create a buffer to receive parent x,y positions
					MPI_Recv(buffer, 3, MPI_INT, 0, mprank+(i*mpisize), MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive parent x,y positions
					localFracElems[i][0] = buffer[1] + displacementX;
					localFracElems[i][1] = buffer[2] + currFracLevel;
					MPI_Send(localFracElems[i], 2, MPI_INT, 0, mprank+(i*mpisize), MPI_COMM_WORLD);
				}
				else if (first)
				{
					if(mpirank == 1)
					{
						localFracElems[i][0] = buffer[1] - displacementX;
						localFracElems[i][1] = buffer[2] + currFracLevel;
						MPI_Send(localFracElems[i], 2, MPI_INT, 0, mprank+(i*mpisize), MPI_COMM_WORLD);
					}
					else
					{
						localFracElems[i][0] = buffer[1] + displacementX;
						localFracElems[i][1] = buffer[2] + currFracLevel;
						MPI_Send(localFracElems[i], 2, MPI_INT, 0, mprank+(i*mpisize), MPI_COMM_WORLD);
					}
					first = false;
				}
				else
				{
					localFracElems[i][0] = buffer[1] - displacementX;
					localFracElems[i][1] = buffer[2] + currFracLevel;
					MPI_Send(localFracElems[i], 2, MPI_INT, 0, mprank+(i*mpisize), MPI_COMM_WORLD);
				}
			}
			maxFracElems = mpthreads; // update max threads for next level control
			currFracLevel++;
		}
	}
  	*/
  }
  printf("Rank %d is ending. \n",mpirank);	
  MPI_Finalize();
  return 0;
}
