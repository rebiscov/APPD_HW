/** Computing the Minimum Spanning Tree of a graph
 * @param N the number of vertices in the graph
 * @param M the number of edges in the graph
 * @param adj the adjacency matrix
 * @param algoName the name of the algorithm to be executed
 */

struct edge{
  int u;
  int v;
};

struct in_tree{
  int u;
  int w;
};

void swap(int *a, int *b){
  int temp;
  temp = *a;
  *a = *b;
  *b = temp;
}

int lexico(int i, int j, int k, int l){ // return 1 if the couple (i, j) < (k, l) for the lexico order (supposing i <= j)
  int temp;
  if (i > j)
    swap(&i, &j);

  if (k > l)
    swap(&k, &l);

  if (i < k || (i == k && j < l))
    return 1;

  return 0;
}


void computeMST(
    int N,
    int M,
    int *adj,
    char *algoName)
{
  int procRank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if (strcmp(algoName, "prim-seq") == 0) { // Sequential Prim's algorithm
    if (procRank == 0) {
      if (numProcs != 1) {
        printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", numProcs);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    // BEGIN IMPLEMENTATION HERE

    int i, j;
    int vertice_set[N];
    struct edge tree[N-1];
    struct in_tree D[N];
    int count = 1; // For now we suppose there is only the vertice 0 in the tree

    memset(vertice_set, 0, sizeof(int)*N);
    vertice_set[0] = 1; // We put 0 in the set of vertice

    for (i = 0; i < N; i++){
      if (adj[i] > 0){
	D[i].u = 0;
	D[i].w = adj[i];
      }
      else
	D[i].w = 0;
    }

    while (count < N){
      struct in_tree candidate;
      int last_i;
      candidate.w = 0;
      for (i = 0; i < N; i++){
	if (!(vertice_set[i]) && D[i].w > 0 && candidate.w == 0){
	  candidate.w = D[i].w;
	  candidate.u = D[i].u;
	  last_i = i;
	}
	else if (!(vertice_set[i]) && D[i].w > 0 && candidate.w > D[i].w) {
	  candidate.w = D[i].w;
	  candidate.u = D[i].u;
	  last_i = i;	  
	}
	else if (!(vertice_set[i]) && D[i].w > 0 && candidate.w == D[i].w) {
	  if (lexico(D[i].u, i, candidate.u, last_i)){
	      candidate.w = D[i].w;
	      candidate.u = D[i].u;
	      last_i = i;	      
	    }
	}
      }
      
      if (last_i < candidate.u)
	printf("%d %d\n", last_i, candidate.u);
      else
	printf("%d %d\n", candidate.u, last_i);

      tree[count - 1].u = last_i;
      tree[count - 1].v = candidate.u;
      vertice_set[last_i] =  1;
      count++;

      for (i = 0; i < N; i++){
	if (!vertice_set[i] && D[i].w > adj[last_i*N + i]){
	  D[i].w = adj[last_i*N + i];
	  D[i].u = last_i;
	}
      }
    }
    

  } else if (strcmp(algoName, "kruskal-seq") == 0) { // Sequential Kruskal's algorithm
    if (procRank == 0) {
      if (numProcs != 1) {
        printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", numProcs);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    // BEGIN IMPLEMENTATION HERE

  } else if (strcmp(algoName, "prim-par") == 0) { // Parallel Prim's algorithm
    // BEGIN IMPLEMENTATION HERE

  } else if (strcmp(algoName, "kruskal-par") == 0) { // Parallel Kruskal's algorithm
    // BEGIN IMPLEMENTATION HERE

  } else { // Invalid algorithm name
    if (procRank == 0) {
      printf("ERROR: Invalid algorithm name: %s.\n", algoName);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}
