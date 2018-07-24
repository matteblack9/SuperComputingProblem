#include "matrixconstructor.h"

void construct_poissonmatrix(int firstrowindex, int firstgrid_x, int lastgrid_x, int numgrid_x, int numgrid_y, double *poissonmatrix)
{
    int count, ii, jj, kk, currentrowindex, localindex;
    int columnsize = numgrid_x*numgrid_y;

    count = 0;

    for(ii=firstgrid_x; ii<=lastgrid_x; ii++)
    {
	for(jj=1; jj<=numgrid_y; jj++)
	{
	    currentrowindex = (ii-1)*numgrid_y + jj-1;
	    localindex = currentrowindex-firstrowindex+1;

	    if(ii==1)
	    {
		if(jj==1)
		{ 
		    poissonmatrix[localindex*columnsize+currentrowindex] = -6.0;
		    poissonmatrix[localindex*columnsize+currentrowindex+1] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex+numgrid_y] = 1.0;
		}
		else if(jj==numgrid_y)
		{
		    poissonmatrix[localindex*columnsize+currentrowindex-1] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex] = -6.0;
		    poissonmatrix[localindex*columnsize+currentrowindex+numgrid_y] = 1.0;
		}
		else
		{
		    poissonmatrix[localindex*columnsize+currentrowindex-1] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex] = -5.0;
		    poissonmatrix[localindex*columnsize+currentrowindex+1] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex+numgrid_y] = 1.0;      
		}
	    }
	    else if(ii==numgrid_x)
	    {
		if(jj==1)
		{
		    poissonmatrix[localindex*columnsize+currentrowindex-numgrid_y] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex] = -6.0;
		    poissonmatrix[localindex*columnsize+currentrowindex+1] = 1.0;
		}
		else if(jj == numgrid_y)
		{
		    poissonmatrix[localindex*columnsize+currentrowindex-numgrid_y] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex-1] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex] = -6.0; 
		}
		else
		{
		    poissonmatrix[localindex*columnsize+currentrowindex-numgrid_y] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex-1] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex] = -5.0;
		    poissonmatrix[localindex*columnsize+currentrowindex+1] = 1.0; 
		}
	    }
	    else{
		if(jj==1)
		{
		    poissonmatrix[localindex*columnsize+currentrowindex-numgrid_y] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex] = -5.0;
		    poissonmatrix[localindex*columnsize+currentrowindex+1] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex+numgrid_y] = 1.0;
		}

		else if(jj==numgrid_y)
		{
		    poissonmatrix[localindex*columnsize+currentrowindex-numgrid_y] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex-1] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex] = -5.0;
		    poissonmatrix[localindex*columnsize+currentrowindex+numgrid_y] = 1.0; 
		}
		else
		{
		    poissonmatrix[localindex*columnsize+currentrowindex-numgrid_y] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex-1] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex] = -4.0;
		    poissonmatrix[localindex*columnsize+currentrowindex+1] = 1.0;
		    poissonmatrix[localindex*columnsize+currentrowindex+numgrid_y] = 1.0; 
		} 
	    }
	}      
    }
}
