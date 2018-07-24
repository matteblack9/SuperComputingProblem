#include<stdio.h>
#include<stdlib.h>
#include<math.h>
void advance_field(int nx, int ny, double *f, double *g) {
	int i, j, idx;
	for (i=1; i<nx-1; i++) {
		for (j=1; j<ny-1; j++) {
			idx = i*ny + j;

			f[idx] = 0.49*(g[idx+ny] + g[idx-ny] + g[idx+1] + g[idx-1] - 4*g[idx]) 
			         + 2*g[idx] - f[idx];
		}
	}
}




int main() {
	const int nx=2000, ny=2000;
	int tmax=500;
	double *f, *g;
	int i, j, idx, tstep;


    //-------------------------------------------------------------------------
    // initialize coefficient and fields
    //-------------------------------------------------------------------------
	f = (double*)malloc(nx*ny*sizeof(double));
	g = (double*)malloc(nx*ny*sizeof(double));

	for (i=0; i<nx; i++) {
		for (j=0; j<ny; j++) { 
			f[i*ny + j] = 0;
			g[i*ny + j] = 0;
		}
	}
	

    //-------------------------------------------------------------------------
  	// main loop for the time evolution
    //-------------------------------------------------------------------------
	for (tstep=1; tstep<=tmax; tstep++) {
    	//------------------------------------
    	// point source
    	//------------------------------------
		idx = (nx/2-1)*ny + (ny/2-1);
	    g[idx] = sin(0.02*tstep);


		advance_field(nx, ny, f, g);
		advance_field(nx, ny, g, f);
	}


    //-------------------------------------------------------------------------
  	// gather fields and save as binary files
    //-------------------------------------------------------------------------
	FILE *fout;
	fout = fopen("field.bin", "wb");
	fwrite(f, sizeof(double), nx*ny, fout);
	fclose(fout);

	return 0;
}
