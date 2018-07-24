#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <inttypes.h>
#include <stdbool.h>

#define BUF_SIZE 2048
#define CHECK_MPI_REQ check_mpi_reqs(startvtx, thislevel, myrank,\
                            nproc, &recvreq_active, &recvreq, recvbuf,\
                            sendreqs_active, sendreqs,\
                            &num_ranks_done, queue, &back,\
                            level, levelsize);

static float value;

typedef struct graphstruct { // A graph in compressed-adjacency-list (CSR) form
	int nv;            // number of vertices
	int64_t ne;            // number of edges
	int *nbr;          // array of neighbors of all vertices
	int *firstnbr;     // index in nbr[] of first neighbor of each vtx
} graph;

int read_graph_from_file(graph **G, int *startvtx, const char* filename);
void print_CSR_graph (const graph *G);
void sendrecv_graph (graph *G, int nproc, int myrank, int **para_range_istart, int **para_range_iend, int *);
void pbfs (const int s, const graph *G, const int nproc, const int myrank, 
		 int *para_range_istart, int *para_range_iend, int *nlevelsp, 
         int **levelsizep);
int para_range(int n1, int n2, int nproc, int myrank, int *istart, int *iend);
int get_ranknum(const int nproc, const int vertex, const int *para_range_istart, const int *para_range_iend, int *rank);

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

void sendrecv_graph (graph *G, int nproc, int myrank, int **para_range_istart, int **para_range_iend, int *startvtx) {
	int i;
	int *sendcounts1, *displs1, *recvbuf1;
	int *sendcounts2, *displs2, *recvbuf2;
	int total_nv, tmp1, istart, iend;
	
	*para_range_istart = (int *)malloc(nproc*sizeof(int));
	*para_range_iend = (int *)malloc(nproc*sizeof(int));
	if (myrank == 0) {
		total_nv = G->nv;
	}
	MPI_Bcast ( startvtx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast ( &total_nv, 1, MPI_INT, 0, MPI_COMM_WORLD);
    value = (float)total_nv/nproc;

	if (myrank == 0) {
		sendcounts1 = (int *)malloc(sizeof(int)*nproc);
		displs1 = (int *)malloc(sizeof(int)*nproc);
		sendcounts2 = (int *)malloc(sizeof(int)*nproc);
		displs2 = (int *)malloc(sizeof(int)*nproc);
	} 
	for (i=0;i<nproc;i++) {	
		para_range(0, total_nv, nproc, i, &istart, &iend);
		(*para_range_istart)[i] = istart;
		(*para_range_iend)[i] = iend;
		if (i!=nproc-1) iend ++;
		tmp1 = iend - istart + 1; 
		if (myrank == 0) {
       		displs1[i] = istart;
        	sendcounts1[i] = tmp1;
			displs2[i] = G->firstnbr[istart];
			sendcounts2[i] = G->firstnbr[iend] - G->firstnbr[istart];
		}
        if (i==myrank) {
			recvbuf1 = (int *) malloc (sizeof(int)*tmp1);
			G->nv = tmp1;
		}
	}
	MPI_Scatterv( G->firstnbr, sendcounts1, displs1, MPI_INT, recvbuf1, G->nv, MPI_INT, 0, MPI_COMM_WORLD); 
	G->ne = recvbuf1[G->nv-1]-recvbuf1[0];
	recvbuf2 = (int *) malloc (sizeof(int)*G->ne);
	MPI_Scatterv( G->nbr, sendcounts2, displs2, MPI_INT, recvbuf2, G->ne, MPI_INT, 0, MPI_COMM_WORLD); 
	if (myrank == 0) {
		free(G->firstnbr);
		free(G->nbr);
		free(sendcounts1);
		free(displs1);
		free(sendcounts2);
		free(displs2);
	}
	G->firstnbr = recvbuf1;
	G->nbr = recvbuf2;
}

static inline void check_mpi_reqs(
				const int startvtx, const int thislevel, const int myrank, 
				const int nproc, bool *recvreq_active, MPI_Request *recvreq, 
				int *recvbuf, bool *sendreqs_active, MPI_Request *sendreqs,
				int *num_ranks_done, int *queue, int *back, 
				int *level, int *levelsize) {
	int i, flag, count;
	while (*recvreq_active) {
		MPI_Status st;
		MPI_Test(recvreq, &flag, &st);
		if (flag) {
			*recvreq_active = false;
			MPI_Get_count(&st, MPI_INT, &count);
			if (count == 0) {
				(*num_ranks_done)++;
			} else {
				for (i=0;i<count;i++) {
					const int w = recvbuf[i];
					if (level[w-startvtx] == -1) { // w has not already been reached
						level[w-startvtx] = thislevel+1; 
						levelsize[thislevel+1]++;
						queue[(*back)++] = w;    // put w on queue to explore
					}
				}
			}
			if (*num_ranks_done < nproc) {
				MPI_Irecv(recvbuf, BUF_SIZE, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, recvreq);
				*recvreq_active=true;
			} 
		} else break;
	}
	for (i=0;i<nproc;i++) {
		if (sendreqs_active[i]) {
			MPI_Test(&sendreqs[i], &flag, MPI_STATUS_IGNORE);
			if (flag) sendreqs_active[i] = false;
		}
	}
}

void pbfs (const int s, const graph *G, const int nproc, const int myrank, 
		 int *para_range_istart, int *para_range_iend, int *nlevelsp, 
         int **levelsizep)  {
	int *level, *levelsize;
	int thislevel;
	int *queue;
	int back, front;
	int i, dest; 
	int startvtx, endvtx;
	bool recvreq_active = false;
	bool *sendreqs_active;
	MPI_Request recvreq;
	MPI_Request *sendreqs;
	int num_ranks_done;
	int *sendbuf, *sendcnt, *recvbuf;

	level = (int *) malloc(G->nv * sizeof(int));
	levelsize = *levelsizep = (int *) calloc(G->nv, sizeof(int));
	queue = (int *) malloc(G->nv * sizeof(int));
	sendcnt = (int *) malloc(nproc * sizeof(int));
	sendreqs = (MPI_Request *) malloc(nproc * sizeof(MPI_Request));
    sendbuf = (int *) malloc(nproc * BUF_SIZE *sizeof(int));
    recvbuf = (int *) malloc(BUF_SIZE *sizeof(int));
	sendreqs_active = malloc(sizeof(bool)*nproc);
	for (i=0;i<nproc;i++) sendreqs_active[i] = false;

	// get start and end vertex
	startvtx = para_range_istart[myrank];
	endvtx = para_range_iend[myrank];

	// initially, queue is empty, all levels are -1
	back = 0;   // position next element will be added to queue
	front = 0;  // position next element will be removed from queue
	for (i=0;i<G->nv;i++) level[i] = -1;

	// assign the starting vertex level 0 and put it on the queue to explore
	thislevel = 0;
	get_ranknum(nproc, s, para_range_istart, para_range_iend, &i);
	if (myrank == i) {
		queue[back++] = s;
		levelsize[0] = 1;
		level[s] = 0;
	}
 
	// loop over levels, then over vertices at this level, then over neighbors
	while(1) {
		levelsize[thislevel+1] = 0;
		memset(sendcnt, 0, nproc*sizeof(int));
		num_ranks_done = 1;

		if (num_ranks_done < nproc) {
			MPI_Irecv(recvbuf, BUF_SIZE, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &recvreq);
			recvreq_active = true;
		}

		for (i = 0; i < levelsize[thislevel]; i++) {
			int v,e;
			int last_w = -1;		
			CHECK_MPI_REQ;
			v = queue[front++];	// v is the current vertex to explore from
			for (e = G->firstnbr[v-startvtx]; e < G->firstnbr[v-startvtx+1]; e++) {
				const int w = G->nbr[e-G->firstnbr[0]];          // w is the current neighbor of v
				if (w == v || last_w==w) continue;
				last_w = w;
				// w exists in current process
				if (w >= startvtx && w <= endvtx) {
               		if (level[w-startvtx] == -1) { // w has not already been reached
                   		level[w-startvtx] = thislevel+1;
                   		levelsize[thislevel+1]++;
                   		queue[back++] = w;    // put w on queue to explore
               		}
				} else { // w exists in remote process
					// get the rank # of w 
					get_ranknum(nproc, w, para_range_istart, para_range_iend, &dest);

					while(sendreqs_active[dest]) 			
						CHECK_MPI_REQ;
					// add (v,w) pair to sendbuf list
                    sendbuf[dest*BUF_SIZE + sendcnt[dest]++] = w;
					if (sendcnt[dest] == BUF_SIZE) {
						MPI_Isend(&sendbuf[dest*BUF_SIZE], BUF_SIZE, MPI_INT, dest, 0, MPI_COMM_WORLD, &sendreqs[dest]);
						sendreqs_active[dest] = true;
						sendcnt[dest] = 0;
					}
				}
			}
		}
		for (i=1;i<nproc;i++) {
			dest = (myrank + i)&(nproc-1);
			if (sendcnt[dest] != 0) {
				while(sendreqs_active[dest]) 			
					CHECK_MPI_REQ;
				MPI_Isend(&sendbuf[dest*BUF_SIZE], sendcnt[dest], MPI_INT, 
							dest, 0, MPI_COMM_WORLD, &sendreqs[dest]);
				sendreqs_active[dest] = true;
				sendcnt[dest] = 0;
			}
			while(sendreqs_active[dest]) 			
				CHECK_MPI_REQ;
			MPI_Isend(&sendbuf[dest*BUF_SIZE], 0, MPI_INT, dest, 0, MPI_COMM_WORLD, &sendreqs[dest]);
			sendreqs_active[dest] = true;
			while(sendreqs_active[dest]) 			
				CHECK_MPI_REQ;
		}
		while (num_ranks_done < nproc) 
			CHECK_MPI_REQ;
		MPI_Allreduce(MPI_IN_PLACE, &levelsize[thislevel],  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		if (levelsize[thislevel] == 0) break;
		thislevel ++;			
	} 
	*nlevelsp = thislevel;

	free(queue);	
	free(sendbuf); 
	free(recvbuf); 
	free(sendcnt); 
	free(sendreqs);
	free(sendreqs_active);

}

int para_range(int n1, int n2, int nproc, int myrank, int *istart, int *iend) {
	int iw1, iw2;
	iw1=(int)(n2-n1+1)/nproc; 
	iw2=(n2-n1+1)%(nproc);
	*istart= myrank*iw1+n1+fmin(myrank,iw2);
	*iend= *istart+iw1-1; 
	if(iw2 > myrank) (*iend)++;
	return 0;
}

int get_ranknum(const int nproc, const int vertex, const int *para_range_istart, const int *para_range_iend, int *rank) {
    int i;
    const int j = (int)floor(vertex/value);
    *rank = -1;
    if (para_range_istart[j] <= vertex && para_range_iend[j] >= vertex) {
        *rank = j;
        return 0;
    } else if (para_range_istart[j] > vertex) {
        for (i=j-1;i>=0;i--) {
            if (para_range_istart[i] <= vertex && para_range_iend[i] >= vertex) {
                *rank = i;
                return 0;
            }
        }
    } else {
        for (i=j+1;i<nproc;i++) {
            if (para_range_istart[i] <= vertex && para_range_iend[i] >= vertex) {
                *rank = i;
                return 0;
            }
        }
    }
    return 1;
}

int main (int argc, char* argv[]) {
	graph *G;
	int *levelsize;
	int nlevels;
	int startvtx=0;
	int i, reached;
	int myrank, nproc;
	const char* filename;
	double starttime, elapsedtime1, elapsedtime2;
	int *para_range_istart, *para_range_iend;

	if (argc == 2) {
		filename = argv[1];
	} else {
		fprintf(stdout, "usage:   pbfs <edgelistfile>\n");
		exit(1);
	}

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	if (myrank == 0) { 
  		starttime = MPI_Wtime();
		if (read_graph_from_file(&G, &startvtx, filename)) return -1;
  		fprintf(stdout, "Elapsed time to read graph data: %f\n", MPI_Wtime() - starttime);
		if (G->nv < nproc) {
			fprintf(stderr, "# of processes should be more than # of vertex:%d\n", G->nv);
            MPI_Abort(MPI_COMM_WORLD, 1);
  			MPI_Finalize();
			exit(-1);
		}
  		print_CSR_graph (G);
	} else {
		G = (graph *) calloc(1, sizeof(graph));
	}

	starttime = MPI_Wtime();
	sendrecv_graph(G, nproc, myrank, &para_range_istart, &para_range_iend, &startvtx);
	MPI_Barrier(MPI_COMM_WORLD);
	elapsedtime1 = MPI_Wtime() - starttime;
	starttime = MPI_Wtime();
	if (myrank == 0) fprintf(stdout, "Starting vertex for BFS is %d.\n\n",startvtx);
	pbfs (startvtx, G, nproc, myrank, para_range_istart, para_range_iend, &nlevels, &levelsize);
	free(para_range_istart);
	free(para_range_iend);
	MPI_Barrier(MPI_COMM_WORLD);
	elapsedtime2 = MPI_Wtime() - starttime;

	if (myrank == 0) {
		reached = 0;
  		for (i = 0; i < nlevels; i++) reached += levelsize[i];
  		fprintf(stdout,"Breadth-first search from vertex %d reached %d levels and %d vertices.\n",
   			startvtx, nlevels, reached);
  		for (i = 0; i < nlevels; i++) fprintf(stdout,"level %d vertices: %d\n", i, levelsize[i]);
    	printf("\n");
		fprintf(stdout, "Elapsed time to scatter graph data: %f\n", elapsedtime1);
    	fprintf(stdout, "Elapsed time to search: %f\n", elapsedtime2);
		fprintf(stdout, "Total elasped time: %f, GTEPs: %f\n", elapsedtime1 + elapsedtime2, (reached/(elapsedtime1 + elapsedtime2)/1000000));
	}

	MPI_Finalize();
	return 0;
}




