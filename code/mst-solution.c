/** Computing the Minimum Spanning Tree of a graph
 * @param N the number of vertices in the graph
 * @param M the number of edges in the graph
 * @param adj the adjacency matrix
 * @param algoName the name of the algorithm to be executed
 */

// structures and functions for the Prim algorithm

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

int lexico(int i, int j, int k, int l){ // return 1 if the couple (i, j) < (k, l) for the lexico order
  int temp;
  if (i > j)
    swap(&i, &j);

  if (k > l)
    swap(&k, &l);

  if (i < k || (i == k && j < l))
    return 1;

  return 0;
}

// structures and functions for the Kruskal algorithm

struct element{ /* This structure will be convinient for the union find */
  int x;
  int nb;
};

struct w_edge{
  int u; /* u and v are vertices, w is the weight of the edge */
  int v;
  int w;
};

int cmp(const void *e1, const void *e2){
  if( (*(struct w_edge*)e1).w - (*(struct w_edge*)e2).w != 0)
    return (*(struct w_edge*)e1).w - (*(struct w_edge*)e2).w;
  else if ( (*(struct w_edge*)e1).u - (*(struct w_edge*)e2).u != 0)
    return (*(struct w_edge*)e1).u - (*(struct w_edge*)e2).u;
  else
    return (*(struct w_edge*)e1).v - (*(struct w_edge*)e2).v;
}

int find(int x, struct element *P){ /* We find, and compress the structure */
  int y, z;

  y = x;

  while (P[y].x != y)
    y = P[y].x;

  while (x != P[x].x){
    z = P[x].x;
    P[x].x = y;
    x = z;
  }

  return y;
}

void unify(int a, int b, struct element *P){ /* We unify (a & b must be representants of their partitions) */
  if (a == b)
    printf("x and y are already in the same partition\n");
  else if (P[a].nb > P[b].nb){
    P[b].x = a;
    P[a].nb += P[b].nb;
  }
  else{
    P[a].x = b;
    P[b].nb += P[a].nb;
  }
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
	if (!vertice_set[i] && adj[last_i*N + i] > 0 && (D[i].w > adj[last_i*N + i] || D[i].w == 0)){
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

    struct w_edge edges[M];
    struct edge tree[N-1];    
    struct element partition[N];
    int i, j, count = 0;

    for (i = 0; i < N; i++)
      for (j = 0; j <= i; j++)
	if (adj[i*N + j] > 0){
	  edges[count].u = j;
	  edges[count].v = i;
	  edges[count++].w = adj[i*N + j];
	}

    qsort(edges, M, sizeof(struct w_edge), cmp);

    for (i = 0; i < N; i++){
      partition[i].x = i;
      partition[i].nb = 1;
    }

    count = 0;
    i = 0;

    while (count < N-1){
      int a = find(edges[i].u, partition), b = find(edges[i].v, partition);
      if (a != b && edges[i].u != edges[i].v){
	tree[count].u = edges[i].u;
	tree[count].v = edges[i].v;
	unify(a, b, partition);
	count++;
      }
      i++;
    }

    for (i = 0; i < N-1; i++)
      printf("%d %d\n", tree[i].u, tree[i].v);

  } else if (strcmp(algoName, "prim-par") == 0) { // Parallel Prim's algorithm
    // BEGIN IMPLEMENTATION HERE

    int vertice_set[N];
    int d_size = procRank != numProcs - 1 ? ceil((float)N/numProcs) : N - ceil((float)N/numProcs)*(numProcs-1);
    int offset = procRank * ceil((float)N / numProcs);
    struct in_tree D[d_size];
    struct edge tree[N-1];
    int i;
    int recv[3*numProcs];
    int send[3];

    memset(vertice_set, 0, sizeof(int)*N);
    vertice_set[0] = 1;

    for (i = 0; i < d_size; i++)
      if (adj[i*N] > 0){
	D[i].w = adj[i*N];
	D[i].u = 0;
      }
      else{
	D[i].w = 0;
	D[i].u = 0;
      }

    for (int count = 0; count < N - 1; count++){
      struct in_tree candidate;
      int last_i;
      candidate.w = 0;
      for (i = 0; i < d_size; i++){
	if (!(vertice_set[offset + i]) && D[i].w > 0 && candidate.w == 0){
	  candidate.w = D[i].w;
	  candidate.u = D[i].u;
	  last_i = i + offset;
	}
	else if (!(vertice_set[offset + i]) && D[i].w > 0 && candidate.w > D[i].w) {
	  candidate.w = D[i].w;
	  candidate.u = D[i].u;
	  last_i = i + offset;	  
	}
	else if (!(vertice_set[offset + i]) && D[i].w > 0 && candidate.w == D[i].w) {
	  if (lexico(D[i].u, i + offset, candidate.u, last_i)){
	    candidate.w = D[i].w;
	    candidate.u = D[i].u;
 	    last_i = i + offset;	      
	  }
	}
      }
      send[0] = last_i;
      send[1] = candidate.u;
      send[2] = candidate.w;

      MPI_Gather(send, 3, MPI_INT, recv, 3, MPI_INT, 0, MPI_COMM_WORLD);

      if (procRank == 0){
	int u = recv[0], v = recv[1], w = recv[2];
	for (int j = 1; j < numProcs; j++)
	  if (w == 0 || (recv[3*j+2] != 0 && recv[3*j + 2] < w || (recv[3*j + 2] == w && lexico(recv[3*j], recv[3*j + 1], u, v)))){
	    u = recv[3*j];
	    v = recv[3*j + 1];
	    w = recv[3*j + 2];
	  }
	send[0] = u;
	send[1] = v;
	send[2] = w;
      }


      MPI_Bcast(&send[0], 3, MPI_INT, 0, MPI_COMM_WORLD);

      if (send[0] > send[1])
	swap(&send[0], &send[1]);
      tree[count].u = send[0];
      tree[count].v = send[1];
      
      int u;
      if (vertice_set[send[0]] == 0){
	u = send[0];
	vertice_set[u] = 1;
      }
      else{
	u = send[1];
	vertice_set[u] = 1;	
      }
      
      for (i = 0; i < d_size; i++){
	if (D[i].w > adj[i*N + u] || D[i].w == 0){
	  D[i].u = u;
	  D[i].w = adj[i*N + u];
	}
      }
    }

    if (procRank == 0){
      for (int count = 0; count < N - 1; count++)
	printf("%d %d\n", tree[count].u, tree[count].v);
    }


  } else if (strcmp(algoName, "kruskal-par") == 0) { // Parallel Kruskal's algorithm
    // BEGIN IMPLEMENTATION HERE

  } else { // Invalid algorithm name
    if (procRank == 0) {
      printf("ERROR: Invalid algorithm name: %s.\n", algoName);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}
