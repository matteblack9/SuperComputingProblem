#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <fcntl.h>
#include <errno.h>
#include <inttypes.h>

typedef struct graphstruct { // A graph in compressed-adjacency-list (CSR) form
	int nv;            // number of vertices
	int64_t ne;            // number of edges
	int *nbr;          // array of neighbors of all vertices
	int *firstnbr;     // index in nbr[] of first neighbor of each vtx
} graph;

int read_graph_from_file(graph **G, int *startvtx, const char* filename);
void print_CSR_graph (const graph *G);
void bfs    (const int s, const graph *G, int *nlevelsp, int **levelsizep);

static inline double get_time() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + tv.tv_usec * 1.e-6;
}

int read_graph_from_file(graph **G, int *startvtx, const char* filename) {

	FILE *fp = fopen (filename, "r");
	if ( fp == NULL ) {
		fprintf(stderr, "error to open file %s: %s\n", filename, strerror(errno));
		return -1;
	}
	*G = (graph *) malloc(sizeof(graph));
	fread(startvtx, 1, sizeof(int), fp);
	fread(&((*G)->ne), 1, sizeof(int64_t), fp);
	fread(&((*G)->nv), 1, sizeof(int), fp);
	(*G)->nbr = (int *) malloc((*G)->ne*sizeof(int));
	(*G)->firstnbr = (int *) malloc(((*G)->nv+1)*sizeof(int));
	fread((*G)->nbr, (*G)->ne, sizeof(int), fp);
	fread((*G)->firstnbr, (*G)->nv+1, sizeof(int), fp);

	fclose(fp);
	return 0;
}

void print_CSR_graph (const graph *G) {
	int vlimit = 20;
	int elimit = 50;
	int e,v;
	fprintf(stdout, "\nGraph has %d vertices and %"PRId64" edges.\n",G->nv,G->ne);
	fprintf(stdout, "firstnbr =");
	if (G->nv < vlimit) vlimit = G->nv;
	for (v = 0; v <= vlimit; v++) printf(" %d",G->firstnbr[v]);
	if (G->nv > vlimit) fprintf(stdout, " ...");
	fprintf(stdout, "\n");
	fprintf(stdout, "nbr =");
	if (G->ne < elimit) elimit = G->ne;
	for (e = 0; e < elimit; e++) printf(" %d",G->nbr[e]);
	if (G->ne > elimit) fprintf(stdout, " ...");
	fprintf(stdout, "\n\n");
}

void bfs   (const int s,
            const graph *G,
            int *nlevelsp,
            int **levelsizep) {
	int *levelsize;
	int thislevel; 
	int *queue, back, front; 
	int i, v, w, e; 
	int *level = (int *) malloc(G->nv*sizeof(int)); 
	levelsize = *levelsizep = (int *) malloc(G->nv*sizeof(int)); 
	queue = (int *) malloc(G->nv*sizeof(int));

	// initially, queue is empty, all levels are -1
	back = 0;   // position next element will be added to queue
	front = 0;  // position next element will be removed from queue
	for (v = 0; v < G->nv; v++) level[v] = -1;

	// assign the starting vertex level 0 and put it on the queue to explore
	thislevel = 0;
	level[s] = 0;
	levelsize[0] = 1;
	queue[back++] = s;

	// loop over levels, then over vertices at this level, then over neighbors
	while (levelsize[thislevel] > 0) {
		levelsize[thislevel+1] = 0;
		for (i = 0; i < levelsize[thislevel]; i++) {
		   	v = queue[front++];       // v is the current vertex to explore from
			for (e = G->firstnbr[v]; e < G->firstnbr[v+1]; e++) {
				w = G->nbr[e];          // w is the current neighbor of v
				if (level[w] == -1) {   // w has not already been reached
					level[w] = thislevel+1;
					levelsize[thislevel+1]++;
					queue[back++] = w;    // put w on queue to explore
				}
			}
		}
		thislevel = thislevel+1;
	}
	*nlevelsp = thislevel;
	free(queue);
	free(level);
}

int main (int argc, char* argv[]) {
	graph *G;
	int *levelsize;
	int nlevels;
	int startvtx;
	int i, reached;
	double starttime, elapsedtime;
	const char* filename;

	if (argc == 2) {
		filename = argv[1];
	} else {
		fprintf(stderr, "usage: bfs <filename>\n");
		exit(1);
	}
	starttime = get_time();
	if (read_graph_from_file(&G, &startvtx, filename)) return -1;
	fprintf(stdout, "Elapsed Time to read and construct graph: %f\n", get_time() - starttime);
	print_CSR_graph (G);

	fprintf(stdout, "Starting vertex for BFS is %d.\n\n",startvtx);
	starttime = get_time();
	bfs (startvtx, G, &nlevels, &levelsize);
	elapsedtime = get_time() -starttime;

	reached = 0;
	for (i = 0; i < nlevels; i++) reached += levelsize[i];
	fprintf(stdout, "Breadth-first search from vertex %d reached %d levels and %d vertices.\n",
		startvtx, nlevels, reached);
	for (i = 0; i < nlevels; i++) printf("level %d vertices: %d\n", i, levelsize[i]);
	fprintf(stdout, "\n");
	fprintf(stdout, "GTEPs: %f\n", (reached/elapsedtime)/1000000);
	fprintf(stdout, "Elapsed Time: %f\n", elapsedtime);

	free(levelsize);
	free(G);
}
