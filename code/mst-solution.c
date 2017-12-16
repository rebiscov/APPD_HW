/** Computing the Minimum Spanning Tree of a graph
 * @param N the number of vertices in the graph
 * @param M the number of edges in the graph
 * @param adj the adjacency matrix
 * @param algoName the name of the algorithm to be executed
 */

#include <stddef.h>

void swap(int *a, int *b){
  int temp;
  temp = *a;
  *a = *b;
  *b = temp;
}

struct edge{
  int u;
  int v;
};

// structures and functions for the Prim algorithm

struct neighbor_in_tree{ /* Represent the closest neighbour u in the tree, w being the "distance" to the tree */
  int u;
  int w;
};

int lexico(int i, int j, int k, int l){ // return 1 if the couple (i, j) < (k, l) for the lexicographic order
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
  int x; /* the "parent" of the element, a represent of a class has himself as his parent */
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

/* Functions for the union-find structure */

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

/* Function for kruskal-par */

struct array_of_forest {
  int size;  
  struct w_edge *forest;
};

struct array_of_forest* kruskal(int N, int numberEdges, struct w_edge *edges){
  int i, count = 0, a, b;
  struct w_edge *forest = malloc(sizeof(struct w_edge)*numberEdges);
  struct element partition[N];
  struct array_of_forest *retvalue = malloc(sizeof(struct array_of_forest));

  qsort(edges, numberEdges, sizeof(struct w_edge), cmp);

  for (i = 0; i < N; i++){
    partition[i].x = i;
    partition[i].nb = 1;
  }

  for (i = 0; i < numberEdges; i++){
    a = find(edges[i].u, partition);
    b = find(edges[i].v, partition);

    if (a != b && edges[i].u != edges[i].v){
      forest[count++] = edges[i];
      unify(a, b, partition);
    }
  }
  retvalue->forest = forest;
  retvalue->size = count;

  return retvalue;
}

struct array_of_forest* merge(int N, struct array_of_forest *F1, struct array_of_forest *F2){
  struct element partition[N];
  int i, j, a, b, count = 0, sumForest = F1->size + F2->size, indexF1 = 0, indexF2 = 0;
  struct array_of_forest *forest = malloc(sizeof(struct array_of_forest));
  struct w_edge *forest1 = F1->forest, *forest2 = F2->forest, candidate;

  forest->forest = malloc(sumForest*sizeof(struct w_edge));

  for (i = 0; i < N; i++){
    partition[i].x = i;
    partition[i].nb = 1;
  }

  for (i = 0; i < sumForest; i++){
    if (indexF2 >= F2->size || cmp(&forest1[indexF1], &forest2[indexF2]) < 0)
      candidate = forest1[indexF1++];
    
    else if (indexF1 >= F1->size || cmp(&forest1[indexF1], &forest2[indexF2]) > 0)
      candidate = forest2[indexF2++];

    else{
      candidate = forest1[indexF1];
      indexF1++;
      indexF2++;
    }
    
    a = find(candidate.u, partition);
    b = find(candidate.v, partition);
    if (a != b){
      (forest->forest)[count++] = candidate;
      unify(a, b, partition);
    }
  }

  forest->size = count;

  return forest;
}

/* Define custom MPI_Op for finding min in prim-par using allreduce  */

void minimumEdge(struct w_edge *in, struct w_edge *inout, int *len, MPI_Datatype *dptr){
  int i;

  for (i = 0; i < *len; i++){
    if (in[i].w == 0 || inout[i].w == 0)
      inout[i] = in[i].w == 0 ? inout[i] : in[i];
    
    else if (in[i].w != inout[i].w)
      inout[i] = in[i].w < inout[i].w ? in[i] : inout[i]; /* put in the result the minimum edge */
    
    else if (lexico(in[i].u, in[i].v, inout[i].u, inout[i].v))
      inout[i] = in[i];
    /* else, do nothing */
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

  MPI_Datatype MPI_WEdge; /* I define a custom type for weighted edge */

  MPI_Type_contiguous(3, MPI_INT, &MPI_WEdge);
  MPI_Type_commit(&MPI_WEdge);

  MPI_Op MPI_Min;
  MPI_Op_create((MPI_User_function*)minimumEdge, 1, &MPI_Min);

  if (strcmp(algoName, "prim-seq") == 0) { // Sequential Prim's algorithm
    if (procRank == 0) {
      if (numProcs != 1) {
        printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", numProcs);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    // BEGIN IMPLEMENTATION HERE

    int i, j;
    int verticeSet[N];
    struct edge tree[N-1];
    struct neighbor_in_tree D[N];
    int count = 1; // For now we suppose there is only the vertice 0 in the tree

    memset(verticeSet, 0, sizeof(int)*N);
    verticeSet[0] = 1; // We put 0 in the set of vertice

    for (i = 0; i < N; i++){
      if (adj[i] > 0){
	D[i].u = 0;
	D[i].w = adj[i];
      }
      else
	D[i].w = 0;
    }

    while (count < N){
      struct w_edge candidate;
      candidate.w = 0;
      for (i = 0; i < N; i++){
	if (!(verticeSet[i]) && D[i].w > 0 && (candidate.w == 0 || candidate.w > D[i].w)){
	  candidate.w = D[i].w;
	  candidate.u = D[i].u;
	  candidate.v = i;
	}
	else if (!(verticeSet[i]) && D[i].w > 0 && candidate.w == D[i].w) {
	  if (lexico(D[i].u, i, candidate.u, candidate.v)){
	      candidate.w = D[i].w;
	      candidate.u = D[i].u;
	      candidate.v = i;	      
	    }
	}
      }
      
      if (candidate.v < candidate.u)
	printf("%d %d\n", candidate.v, candidate.u);
      else
	printf("%d %d\n", candidate.u, candidate.v);

      tree[count - 1].u = candidate.v;
      tree[count - 1].v = candidate.u;
      verticeSet[candidate.v] =  1;
      count++;

      for (i = 0; i < N; i++){
	if (!verticeSet[i] && adj[candidate.v*N + i] > 0 && (D[i].w > adj[candidate.v*N + i] || D[i].w == 0)){
	  D[i].w = adj[candidate.v*N + i];
	  D[i].u = candidate.v;
	}
      }
    }

    /* END OF PRIM-SEQ */

  } else if (strcmp(algoName, "kruskal-seq") == 0) { // Sequential Kruskal's algorithm
    if (procRank == 0) {
      if (numProcs != 1) {
        printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", numProcs);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    // BEGIN IMPLEMENTATION HERE

    struct w_edge edges[M]; /* the set of edges that will be sorted */
    struct edge tree[N-1]; /* the final tree */
    struct element partition[N]; /* structure to keep track of the partitions, we do not want to create cycle */
    int i, j, count = 0;

    for (i = 0; i < N; i++) /* We add all the edges to the array edges */
      for (j = 0; j <= i; j++)
	if (adj[i*N + j] > 0){
	  edges[count].u = j;
	  edges[count].v = i;
	  edges[count++].w = adj[i*N + j];
	}

    qsort(edges, M, sizeof(struct w_edge), cmp); /* we sort them */

    for (i = 0; i < N; i++){ /* intialization of the union-find structure, at the beginnig, all the vertices are alone in their partitions */
      partition[i].x = i; /* partition[i] points to the representant of i, for now it is i  */
      partition[i].nb = 1; /* i is alone in its partition */
    }

    count = 0;
    i = 0;

    while (count < N-1){ /* while we do not have a tree... */
      int a = find(edges[i].u, partition), b = find(edges[i].v, partition); /* check if we can add the edge (u, v) without creating a cycle */
      
      if (a != b && edges[i].u != edges[i].v){
	tree[count].u = edges[i].u;
	tree[count].v = edges[i].v;
	unify(a, b, partition); /* put a and b in the same set */
	count++;
      }
      i++;
    }

    for (i = 0; i < N-1; i++)
      printf("%d %d\n", tree[i].u, tree[i].v);

    /* END OF KRUSKAL-SEQ */

  } else if (strcmp(algoName, "prim-par") == 0) { // Parallel Prim's algorithm
    // BEGIN IMPLEMENTATION HERE

    int verticeSet[N]; /* boolean array to keep track of the vertices in the array */
    int numRows = procRank != numProcs - 1 ? ceil((float)N/numProcs) : N - ceil((float)N/numProcs)*(numProcs-1); /* Number of rows that the processor has */
    int offset = procRank * ceil((float)N / numProcs);
    struct neighbor_in_tree D[numRows];
    struct edge tree[N-1];
    int i;

    memset(verticeSet, 0, sizeof(int)*N);
    verticeSet[0] = 1;

    for (i = 0; i < numRows; i++) /* I initialize the D vector, in addition of the minimun distance, I add the vertice for which this distance is tight */
      if (adj[i*N + 0] > 0){
	D[i].w = adj[i*N];
	D[i].u = 0;
      }
      else{
	D[i].w = 0;
	D[i].u = 0;
      }

    for (int count = 0; count < N - 1; count++){
      struct w_edge candidate; /* I store the candidate for the edge that minimizes the distance from a vertex to the tree */
      candidate.w = 0;
      for (i = 0; i < numRows; i++){
	if (!(verticeSet[offset + i]) && D[i].w > 0 && (candidate.w == 0 || candidate.w > D[i].w)){
	  candidate.w = D[i].w;
	  candidate.u = D[i].u;
	  candidate.v = i + offset;
	}
	else if (!(verticeSet[offset + i]) && D[i].w > 0 && candidate.w == D[i].w) {
	  if (lexico(D[i].u, i + offset, candidate.u, candidate.v)){
	    candidate.w = D[i].w;
	    candidate.u = D[i].u;
 	    candidate.v = i + offset;	      
	  }
	}
      }
      struct w_edge recv;
      
      MPI_Allreduce(&candidate, &recv, 1, MPI_WEdge, MPI_Min, MPI_COMM_WORLD); /* We take the edge with minimum weight and we send it to everybody  */

      if (recv.u > recv.v) /* We add the chosen edge in the tree*/
	swap(&recv.u, &recv.v);
      tree[count].u = recv.u;
      tree[count].v = recv.v;
      
      int u;
      if (verticeSet[recv.u] == 0){ /* We add the new vertex to the tree */
	u = recv.u;
	verticeSet[u] = 1;
      }
      else{
	u = recv.v;
	verticeSet[u] = 1;	
      }
      
      for (i = 0; i < numRows; i++){ /* Each processor updates its array D */
	if ((D[i].w > adj[i*N + u] && adj[i*N + u] > 0) || D[i].w == 0){
	  D[i].u = u;
	  D[i].w = adj[i*N + u];
	}
      }
    }

    if (procRank == 0){ /* Processor 0 displays the result */
      for (int count = 0; count < N - 1; count++)
	printf("%d %d\n", tree[count].u, tree[count].v);
    }

    /* END OF PRIM-PAR */

  } else if (strcmp(algoName, "kruskal-par") == 0) { // Parallel Kruskal's algorithm
    // BEGIN IMPLEMENTATION HERE

    struct w_edge edges[M];
    int numRows = procRank != numProcs - 1 ? ceil((float)N/numProcs) : N - ceil((float)N/numProcs)*(numProcs-1); /* Number of rows that the processor has */    
    int offset = procRank * ceil((float)N / numProcs);
    int count = 0;

    for(int i = 0; i < numRows; i++) /* Put the edges in an array, it is not necessary to take all the edges in the submatrix adj, If we look at it globally, we take half the matrix, and because the matrix is symmetric, all the edges are covered by at least one processor */
      for (int j = i + offset; j < N; j++)
	if (adj[i*N + j] > 0){
	  edges[count].u = i + offset;
	  edges[count].v = j;
	  edges[count++].w = adj[i*N + j];
	}

    struct array_of_forest *forest = kruskal(N, count, edges);

    int aliveProcs = numProcs; /* processors still running */
    int n = 0; /* shift for logical shift */
    int currentRank = procRank; /* rank convinient to choose who has to send/recv */

    struct array_of_forest recv;
    
    while(aliveProcs != 1){
      /* If there is an odd number of alive processors, the last one should not receive nor send informations to others */
      if (!(aliveProcs % 2 == 1 && procRank == aliveProcs - 1)){

	if (currentRank % 2 == 0){ /* even processors receive and merge */
	  int src = procRank + (1 << n);
	  MPI_Recv(&recv.size, 1, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	  MPI_Recv(recv.forest, recv.size, MPI_WEdge, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	  struct w_edge *temp = forest->forest;
	  
	  forest = merge(N, forest, &recv);
	  
	  free(recv.forest);
	  free(temp);
	}
	else{ /* odd processors send all and stop */
	  int dest = procRank - (1 << n);
	  MPI_Send(&forest->size, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
	  MPI_Send(forest->forest, forest->size, MPI_WEdge, dest, 0, MPI_COMM_WORLD);

	  free(forest->forest);
	  break; /* end of the work for this processor */
	}
	
      }

      n++;
      aliveProcs = (aliveProcs >> 1) + (aliveProcs & 1);
      currentRank = (currentRank >> 1) + (currentRank & 1);
    }

    if (procRank == 0)
      for (int i = 0; i < forest->size; i++)
	if (forest->forest[i].u < forest->forest[i].v)
	  printf("%d %d\n", forest->forest[i].u, forest->forest[i].v);
	else
	  printf("%d %d\n", forest->forest[i].v, forest->forest[i].u);
    
  } else { // Invalid algorithm name
    if (procRank == 0) {
      printf("ERROR: Invalid algorithm name: %s.\n", algoName);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}
