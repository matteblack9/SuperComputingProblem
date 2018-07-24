¹®Á¦ 1¹ø (C)

#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<float.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

void axmb(int m, int n, double *uu, double *ww, double *xg, double *yg){
	int i,j;
	double pi, hh, hh2;
	hh = xg[1] - xg[0]; hh2 = hh*hh;
	for(i=0;i<m*n;i++) ww[i] = 0.L;

	for(i=1;i<m-1;i++){
		for(j=0;j<n;j++){
			ww[j+n*i] += (-2.L*uu[j+n*i]+uu[j+n*(i-1)] + uu[j+n*(i+1)])/hh2;
		}
	}
	hh = yg[1] - yg[0]; hh2 = hh*hh;
	for(i=1;i<m-1;i++){
		for(j=1;j<n-1;j++){
			ww[j+n*i] += (-2.L*uu[j+n*i]+uu[j-1+n*i] + uu[j+1+n*i])/hh2;
		}
	}
	pi = M_PI;
	for(i=0;i<m;i++){
		for(j=0;j<n;j++){
			ww[j+n*i] -= (-2.l*pi*pi*cos(pi*xg[i])*cos(pi*yg[j]));
		}
	}
}
void vndb(int m, int n, double *uu){
	int i,j;
	for(j=0;j<n;j++) uu[j] = uu[j+1*n];
	for(j=0;j<n;j++) uu[j+n*(m-1)] = uu[j+n*(m-2)];
	for(i=0;i<m;i++) uu[n*i] = uu[1+n*i];
	for(i=0;i<m;i++) uu[n-1+n*i] = uu[n-2+n*i];
}


int main(int argc, char **argv){
	double *uu, *vv, *ww;
	double *xg, *yg;
	int i,j,k,miter,iter;
	double pi, test0, hh,hh2;

	int m,n;

	n = m = 1000;
	xg = (double*)malloc(sizeof(double)*m);
	yg = (double*)malloc(sizeof(double)*n);
	uu = (double*)malloc(sizeof(double)*m*n);
	vv = (double*)malloc(sizeof(double)*m*n);
	ww = (double*)malloc(sizeof(double)*m*n);


	for(i=0;i<m;i++) xg[i] = 0.L + (1.L-0.L)/(double)(m-1)*(double)(i);
	for(j=0;j<n;j++) yg[j] = 0.L + (1.L-0.L)/(double)(n-1)*(double)(j);


	pi = M_PI;
	hh = xg[1] - xg[0];
	hh2 = hh * hh;
	for(i=0;i<m*n;i++) uu[i] = 0.l;
	vndb(m,n,uu);
	for(i=0;i<m*n;i++) vv[i] = uu[i];
	miter = 1000000000;
	miter = 10000;

	for(iter=1;iter<=miter;iter++){
		vndb(m,n,uu);
		axmb(m,n,uu,ww,xg,yg);
		vndb(m,n,ww);
		for(i=1;i<m-1;i++){
			for(j=1;j<n-1;j++){
				uu[j+n*i] = vv[j+n*i] + ww[j+n*i]*(hh2/4.L)*0.9L;
			}
		}
		test0 = -DBL_MAX;
		double amax = -DBL_MAX;
		for(i=0;i<m*n;i++){
			test0 = max(fabs(vv[i]-uu[i]),test0);
			amax = max(fabs(uu[i])+1.E-8L,amax);
		}
		test0 = test0/amax;
		for(i=0;i<m*n;i++){
			vv[i] = uu[i];
		}
		if((iter%1000) == 1 || iter < 1000) printf("%d %24.15e\n",iter, test0);
		if(test0<1.E-6L ) break;
	}

	printf("%24.15g %24.15g\n",uu[0],uu[m-1]);
	printf("%24.15g %24.15g\n",uu[n-1],uu[n-1+n*(m-1)]);
	FILE *wp = fopen("fort.11", "w");
	for(i=0;i<m;i++){
		for(j=0;j<n;j++){
			fprintf(wp,"%g %g %g\n",xg[i],yg[i],uu[j+n*i]);
		}
		fprintf(wp,"\n");
	}
	fprintf(wp,"\n");
	free(xg);
	free(yg);
	free(uu);
	free(ww);
	free(vv);
}
