#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <math.h>

#include "mst-solution.c"

void readGraph(
    char *fileName,
    int *pnVtx,
    int *pnEdge,
    int **padj)
{
  FILE *file = fopen(fileName, "r");
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (file == NULL) {
    if (rank == 0) {
      printf("ERROR: Unable to open the file %s.\n", fileName);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  int nVtx, nEdge;
  int *adj;
  int *tempadj;
  fscanf(file, " %d %d", pnVtx, pnEdge);
  nVtx = *pnVtx; nEdge = *pnEdge;
  int *edgeLeft = NULL, *edgeRight = NULL, *weightL = NULL;
  if (rank == 0)
  {
      edgeLeft = (int *) malloc(nEdge * sizeof(edgeLeft[0]));
      edgeRight = (int *) malloc(nEdge * sizeof(edgeRight[0]));
      weightL = (int *) malloc(nEdge * sizeof(weightL[0]));
  }
  tempadj = (int *) malloc(nVtx * nVtx * sizeof(tempadj[0]));

  //Allocate the matrix and fill it with input data
  int nb_elements = nVtx*(int)ceil((float)nVtx/(float)size);
  if (rank == size-1)
      *padj = adj = (int *) malloc((nVtx*nVtx-(size-1)*nb_elements)*sizeof(adj[0]));
  else
      *padj = adj = (int *) malloc(nb_elements*sizeof(adj[0]));
  
  if (rank == 0)
  {
      int i,j;
      for (i = 0; i < nVtx; i++)
        for (j = 0; j < nVtx; j++)
            tempadj[i*nVtx+j] = 0;

      for (i = 0; i < nEdge; i++) {
        fscanf(file, " %d %d %d", edgeLeft + i, edgeRight + i, weightL + i);
        tempadj[edgeLeft[i]*nVtx + edgeRight[i]] = weightL[i];
        tempadj[edgeRight[i]*nVtx + edgeLeft[i]] = weightL[i];
      }
      free(edgeLeft);
      free(edgeRight);
      free(weightL);

  }

  if (size > 1)
  {
    //  MPI_Bcast(adj, nVtx*nVtx, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Comm commFirst;
    MPI_Comm_split(MPI_COMM_WORLD,rank/(size-1),rank,&commFirst);

    //SCATTER THE MATRIX FOR P-1 PROCESS
    if (rank != size-1)
      MPI_Scatter(tempadj, nb_elements, MPI_INT, adj, nb_elements, MPI_INT, 0, commFirst);

    //SEND THE LAST ELEMENTS TO THE LAST PROCESS
    if (rank == 0)
      MPI_Send(tempadj+(size-1)*nb_elements,nVtx*nVtx-(size-1)*nb_elements,MPI_INT,size-1,0,MPI_COMM_WORLD);
    else if (rank == size-1)
      MPI_Recv(adj,nVtx*nVtx-(size-1)*nb_elements,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

  } else {
    memcpy(adj,tempadj,nVtx*nVtx*sizeof(tempadj[0]));
  }

  free(tempadj);

  fclose(file);
}

void printUsage() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    printf("Usage: mpirun -np [num-procs] ./mst [graph-file-name] [algo-name]\n"
        "Arguments:\n"
        "\t[num-procs]: The number of MPI ranks to be used. For the sequential algorithm it has to be 1.\n"
        "\t[graph-file-name]: Name of the graph file.\n"
        "\t[bisect-algo-name]: Name of the graph bisection algorithm. There are three possibilities:\n"
        "\t\tprim-seq: Sequential Prim's algorithm.\n"
        "\t\tkruskal-seq: Sequential Kruskal's algorithm.\n"
        "\t\tprim-par: Parallel Prim's algorithm.\n"
        "\t\tkruskal-par: Parallel Kruskal's algorithm.\n"
        );
  }
}


int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  if (argc < 3) {
    printUsage();
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int procRank,commSize;
  MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
  MPI_Comm_size(MPI_COMM_WORLD, &commSize);
  char *algoName = argv[2];
  int nVtx, nEdge;
  int *adj = NULL;
  int seed;
  if (argc == 6) { seed = atoi(argv[5]); }
  else { seed = time(0); }

  readGraph(argv[1], &nVtx, &nEdge, &adj);
  MPI_Barrier(MPI_COMM_WORLD);

  double startTime = MPI_Wtime();
  computeMST(nVtx, nEdge, adj, algoName);

  if (procRank == 0) { printf("computeMST took %e seconds.\n", MPI_Wtime() - startTime); }

  MPI_Finalize();

  free(adj);
  return 0;
}
