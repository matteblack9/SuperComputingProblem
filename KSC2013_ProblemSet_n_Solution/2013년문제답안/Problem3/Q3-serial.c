#include<stdio.h>
#include<stdlib.h>
#include<math.h>



void advance_field(int nx, int ny, double *f, double *g, double *c) {
	int i, j, idx;

	for (i=1; i<nx-1; i++) {
		for (j=1; j<ny-1; j++) {
			idx = i*ny + j;

			f[idx] = c[idx]*(g[idx+ny] + g[idx-ny] + g[idx+1] + g[idx-1] - 4*g[idx]) 
		             + 2*g[idx] - f[idx];
		}
	}
}



void periodic_y(int nx, int ny, double *f) {
	int i, j;

	for (i=0; i<nx; i++) {
		f[i*ny+(ny-1)] = f[i*ny+1];
		f[i*ny+0] = f[i*ny+(ny-2)];
	}
}



int main() {
	int nx=2000, ny=2000;
	int tmax=1500;
	int width=120, thick=60, gap=1200;	// slit parameters
	int distance=800;
	double c0=0.49;
	double *f, *g, *c;
	int i, j, tstep;


    //-------------------------------------------------------------------------
    // initialize coefficient and fields
    //-------------------------------------------------------------------------
	f = (double*)malloc(nx*ny*sizeof(double));
	g = (double*)malloc(nx*ny*sizeof(double));
	c = (double*)malloc(nx*ny*sizeof(double));

	for (i=0; i<nx*ny; i++) {
		f[i] = 0;
		g[i] = 0;
		c[i] = c0;
	}
	

    //-------------------------------------------------------------------------
  	// slit structure
    //-------------------------------------------------------------------------
	for (i=nx/2-1; i<nx/2+thick; i++) {
		for (j=0; j<ny; j++)
			c[i*ny+j] = 0;
	}

	for (i=nx/2-1; i<nx/2+thick; i++) {
		for (j=ny/2-(gap+width)/2-1; j<ny/2+(gap+width)/2; j++)
			c[i*ny+j] = c0;
	}

	for (i=nx/2-1; i<nx/2+thick; i++) {
		for (j=ny/2-(gap-width)/2-1; j<ny/2+(gap-width)/2; j++)
			c[i*ny+j] = 0;
	}


    //-------------------------------------------------------------------------
  	// main loop for the time evolution
    //-------------------------------------------------------------------------
	for (tstep=1; tstep<=tmax; tstep++) {
    	//------------------------------------
    	// point source
    	//------------------------------------
	    //g[nx/2-1][ny/2-1] = sin(0.02*tstep);


    	//------------------------------------
    	// line source
    	//------------------------------------
		for (j=0; j<ny; j++)
			g[(nx/2-distance-1)*ny+j] = sin(0.02*tstep);


		advance_field(nx, ny, f, g, c);
		periodic_y(nx, ny, f);

		advance_field(nx, ny, g, f, c);
		periodic_y(nx, ny, g);
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
