#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<float.h>
#include<mpi.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))



void axmb(int m, int n, int ns1, int ns2, double *uu, double *ww, double *xg, double *yg,
	  int jsta2, int jend1){
  int i,j;
  double pi, hh, hh2;
  hh = xg[1] - xg[0]; hh2 = hh*hh;
  for(i=0;i<m*(ns2-ns1+1);i++) ww[i] = 0.L;
  
  for(j=1;j<ns2-ns1;j++){
    for(i=1;i<m-1;i++){
      ww[i+m*j] += (-2.L*uu[i+m*j]+uu[i-1+m*j] + uu[i+1+m*j])/hh2;
    }
  }
  hh = yg[1] - yg[0]; hh2 = hh*hh;
  for(j=1;j<ns2-ns1;j++){
    for(i=1;i<m-1;i++){
      ww[i+m*j] += (-2.L*uu[i+m*j]+uu[i+m*(j-1)] + uu[i+m*(j+1)])/hh2;
    }
  }
  pi = M_PI;
  for(j=0;j<ns2-ns1+1;j++){
    for(i=0;i<m;i++){
      ww[i+m*j] -= (-2.l*pi*pi*cos(pi*xg[i])*cos(pi*yg[j+ns1-1]));
    }
  }
}
void vnbd(int m, int n, int ns1, int ns2,double *uu){
  int i,j;
  for(j=0;j<ns2-ns1+1;j++) uu[j*m] = uu[1+j*m];
  for(j=0;j<ns2-ns1+1;j++) uu[m-1+m*j] = uu[m-2+m*j];
  if(ns1== 1) for(i=0;i<m;i++) uu[i+m*(1-ns1)] = uu[i+m*(2-ns1)];
  if(ns2 ==n) for(i=0;i<m;i++) uu[i+m*(n-ns1)] = uu[i+m*(n-ns1-1)];
}


void para_range(int n1, int n2, int nid, int myid, int *ista, int *iend){
  int iwork1, iwork2;
  
  iwork1 = (n2-n1+1)/nid;
  iwork2 = (n2-n1+1)%nid;
  *ista = myid*iwork1 + n1 + min(myid, iwork2);
  *iend = *ista + iwork1 - 1;
  if(iwork2 > myid) *iend += 1;
}
void para_range1(int n1, int n2, int nid, int myid, int *ista, int *iend){
  int iwork;
  iwork = (n2-n1)/nid + 1;
  *ista = min(myid * iwork + n1, n2+1);
  *iend = min( *ista + iwork-1, n2);
}


int main(int argc, char **argv){
  double *uu, *vv, *ww;
  double *xg, *yg;
  int i,j,k,miter,iter;
  double pi, test0, hh,hh2;
  int ns1,ns2;
  double *uurow;
  int jsta,jend,jsta2,jend1,iprev,inext;
  MPI_Request isend1,isend2,irecv1,irecv2;
  int iplot,myid,nid;
  MPI_Status istatus;
  double test1,wtime1,wtime2;
  char lexit;
  
  
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nid);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  
  wtime1 = MPI_Wtime();
  
  int m,n;
  
  
  
  n = (m = 2000);
  xg = (double*)malloc(sizeof(double)*m);
  yg = (double*)malloc(sizeof(double)*m);
  
  
  for(i=0;i<m;i++) xg[i] = 0.L + (1.L-0.L)/(double)(m-1)*(double)(i);
  for(j=0;j<n;j++) yg[j] = 0.L + (1.L-0.L)/(double)(n-1)*(double)(j);
  

  pi = M_PI;
  hh = xg[1] - xg[0];
  hh2 = hh * hh;
  
  
  
  para_range(1,n,nid,myid,&jsta,&jend);
  jsta2 = jsta;
  jend1 = jend;
  
  if(myid == 0) jsta2 = 2;
  if(myid == nid-1) jend1 = n-1;
  ns1 = max(1, jsta-1); 
  ns2 = min(n, jend+1);
  
  uu = (double*)malloc(sizeof(double)*m*(ns2-ns1+1));
  vv = (double*)malloc(sizeof(double)*m*(ns2-ns1+1));
  ww = (double*)malloc(sizeof(double)*m*(ns2-ns1+1));
  inext = myid + 1;
  iprev = myid - 1;
  if(myid == nid-1) inext = MPI_PROC_NULL;
  if(myid == 0) iprev = MPI_PROC_NULL;
  
  for(j=jsta; j<=jend;j++) for(i=0;i<m;i++) uu[i+m*(j-ns1)] = 0;
  
  vnbd(m,n,ns1,ns2,uu);
  for(i=0;i<m*(ns2-ns1+1);i++) vv[i] = uu[i];
  miter = 1000000000;
  miter = 10000;
  
  for(iter=1;iter<=miter;iter++){
    vnbd(m,n,ns1,ns2,uu);
    MPI_Isend(uu+m*(jend-ns1), m , MPI_DOUBLE, inext, 1, MPI_COMM_WORLD, &isend1);
    MPI_Isend(uu+m*(jsta-ns1), m , MPI_DOUBLE, iprev, 1, MPI_COMM_WORLD, &isend2);
    MPI_Irecv(uu+m*(jsta-1-ns1), m , MPI_DOUBLE, iprev, 1, MPI_COMM_WORLD, &irecv1);
    MPI_Irecv(uu+m*(jend+1-ns1), m , MPI_DOUBLE, inext, 1, MPI_COMM_WORLD, &irecv2);
    MPI_Wait(&isend1,&istatus);
    MPI_Wait(&isend2,&istatus);
    MPI_Wait(&irecv1,&istatus);
    MPI_Wait(&irecv2,&istatus);
    
    axmb(m,n,ns1,ns2,uu,ww,xg,yg, jsta2, jend1);
    vnbd(m,n,ns1,ns2,ww);
    for(j=jsta2;j<=jend1;j++){
      for(i=0;i<m;i++){
	uu[i+m*(j-ns1)] = vv[i+m*(j-ns1)] + ww[i+m*(j-ns1)]*(hh2/4.L)*0.9L;
      }
    }
    test0 = -1.e30L;
    double amax = -1.e20L,amax1;
    for(i=0;i<m*(ns2-ns1+1);i++){
      test0 = max(fabs(vv[i]-uu[i]),test0);
      amax = max(fabs(uu[i]),amax);
    }
    amax = amax + 1.e-8L;
    MPI_Reduce(&test0, &test1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&amax, &amax1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(myid==0){
      amax = amax1; test0 = test1;
    }
    MPI_Bcast(&amax,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&test0,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
    test0 = test0/amax;
    for(i=0;i<m*(ns2-ns1+1);i++){
      vv[i] = uu[i];
    }
    if(myid ==0){
      if((iter%1000) == 1 || iter < 1000) printf("%d %24.15e\n",iter, test0);
      lexit = 'N';
      if(test0<1.e-6L) lexit = 'Y';
    }
    MPI_Bcast(&lexit, 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    if(lexit == 'Y') break;
  }
  
  if(myid==0) printf("%24.15g %24.15g\n",uu[0 + m *(1-ns1)],uu[m-1 + m *(1-ns1)]);
  if(myid==nid-1) printf("%24.15g %24.15g\n",uu[0 + m * (n-ns1)],uu[m-1 + m *(n-ns1)]);
  /*
    FILE *wp = fopen("sol","w");
    long ioffset = (jsta-1) * sizeof(double)*(m);
    fseek(wp, ioffset, SEEK_SET);
    for(j=jsta;j<=jend;j++){
    fwrite(uu+(j-ns1)*m, sizeof(double), m, wp);
    }
    fclose(wp);
  */
  
  MPI_Barrier(MPI_COMM_WORLD);
  wtime2 = MPI_Wtime();
  if(myid ==0) printf("[C] np=%d \t %g (sec)\n",nid, wtime2-wtime1);
  iplot = 0;
  iplot = 1;
  /*
    if(iplot ==1 && myid == 0) {
    FILE *fp = fopen("sol", "r");
    wp = fopen("fort.11","w");
    
    
    uurow = (double*)malloc(sizeof(double)*m);
    
    for(j=0;j<n;j++){
    fread(uurow,sizeof(double), m,fp);
    for(i=0;i<m;i++){
    fprintf(wp,"%24.15g %24.15g %24.15g\n",xg[i],yg[j],uurow[i]);
    }
    fprintf(wp,"\n");
    }
    fclose(fp);
    fclose(wp);
    free(uurow);
    }
  */
  free(xg);
  free(yg);
  free(uu);
  free(ww);
  free(vv);
  MPI_Finalize();
}
