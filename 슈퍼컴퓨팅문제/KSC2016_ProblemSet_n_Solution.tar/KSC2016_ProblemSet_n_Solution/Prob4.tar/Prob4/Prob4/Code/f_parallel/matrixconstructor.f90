module matrixconstructor
    
contains

    subroutine construct_poissonmatrix(firstrowindex, firstgrid_y, lastgrid_y, numgrid_x, numgrid_y, poissonmatrix)
    
    implicit none

    integer(4) :: firstrowindex, firstgrid_y, lastgrid_y
    integer(4) :: numgrid_x, numgrid_y
    real(8) :: poissonmatrix(numgrid_x*numgrid_y*numgrid_x*numgrid_y)

    integer(4) :: count, ii, jj, kk, currentrowindex, localindex
    integer(4) :: columnsize

    count = 0
    columnsize = numgrid_x * numgrid_y

    do jj = firstgrid_y,lastgrid_y
        do ii = 1, numgrid_x
            currentrowindex = (jj-1)*numgrid_x + ii
            localindex = currentrowindex-firstrowindex
            if(jj.EQ.1) then
                if(ii.EQ.1) then
                    poissonmatrix(localindex*columnsize+currentrowindex) = -6.0
                    poissonmatrix(localindex*columnsize+currentrowindex+1) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex+numgrid_x) = 1.0
                else if(ii.EQ.numgrid_x) then
                    poissonmatrix(localindex*columnsize+currentrowindex-1) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex) = -6.0
                    poissonmatrix(localindex*columnsize+currentrowindex+numgrid_x) = 1.0
                else
                    poissonmatrix(localindex*columnsize+currentrowindex-1) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex) = -5.0
                    poissonmatrix(localindex*columnsize+currentrowindex+1) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex+numgrid_x) = 1.0
                endif
            else if(jj.EQ.numgrid_y) then
                if(ii.EQ.1) then
                    poissonmatrix(localindex*columnsize+currentrowindex-numgrid_x) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex) = -6.0
                    poissonmatrix(localindex*columnsize+currentrowindex+1) = 1.0
                else if(ii.EQ.numgrid_x) then
                    poissonmatrix(localindex*columnsize+currentrowindex-numgrid_x) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex-1) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex) = -6.0
                else
                    poissonmatrix(localindex*columnsize+currentrowindex-numgrid_x) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex-1) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex) = -5.0
                    poissonmatrix(localindex*columnsize+currentrowindex+1) = 1.0
                endif
            else
                if(ii.EQ.1) then
                    poissonmatrix(localindex*columnsize+currentrowindex-numgrid_x) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex) = -5.0
                    poissonmatrix(localindex*columnsize+currentrowindex+1) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex+numgrid_x) = 1.0
                else if (ii.EQ.numgrid_x) then
                    poissonmatrix(localindex*columnsize+currentrowindex-numgrid_x) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex-1) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex) = -5.0
                    poissonmatrix(localindex*columnsize+currentrowindex+numgrid_x) = 1.0
                else
                    poissonmatrix(localindex*columnsize+currentrowindex-numgrid_x) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex-1) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex) = -4.0
                    poissonmatrix(localindex*columnsize+currentrowindex+1) = 1.0
                    poissonmatrix(localindex*columnsize+currentrowindex+numgrid_x) = 1.0
                endif
            endif
        enddo
    enddo

    end subroutine construct_poissonmatrix
    
end module matrixconstructor
    
