subroutine construct_poissonmatrix(firstrowindex, firstgrid_x, lastgrid_x, numgrid_x, numgrid_y, poissonmatrix)

      implicit none

      integer(4) :: firstrowindex, firstgrid_x, lastgrid_x, numgrid_x, numgrid_y
      real(8) :: poissonmatrix( (lastgrid_x - firstgrid_x + 1) * numgrid_y )

      integer(4) :: count, ii, jj, kk, currentrowindex, localindex
      integer(4) :: columnsize

      columnsize = numgrid_x * numgrid_y
      count = 0

      do ii=firstgrid_x, lastgrid_x
          do jj=1, numgrid_y
              currentrowindex = (ii-1)*numgrid_y + jj
              localindex = currentrowindex-firstrowindex

              if(ii.EQ.1) then
                  if(jj.EQ.1) then
                      poissonmatrix(localindex*columnsize+currentrowindex)           = -4.0
                      poissonmatrix(localindex*columnsize+currentrowindex+1)         =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex+numgrid_y) =  1.0
                  else if(jj.EQ.numgrid_y) then
                      poissonmatrix(localindex*columnsize+currentrowindex-1)         =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex)           = -4.0
                      poissonmatrix(localindex*columnsize+currentrowindex+numgrid_y) =  1.0
                  else
                      poissonmatrix(localindex*columnsize+currentrowindex-1)         =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex)           = -4.0
                      poissonmatrix(localindex*columnsize+currentrowindex+1)         =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex+numgrid_y) =  1.0
                  endif 
              else if(ii.EQ.numgrid_x) then
                  if(jj.EQ.1) then
                      poissonmatrix(localindex*columnsize+currentrowindex-numgrid_y) =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex)           = -4.0
                      poissonmatrix(localindex*columnsize+currentrowindex+1)         =  1.0
                  else if(jj.EQ.numgrid_y) then
                      poissonmatrix(localindex*columnsize+currentrowindex-numgrid_y) =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex-1)         =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex)           = -4.0
                  else
                      poissonmatrix(localindex*columnsize+currentrowindex-numgrid_y) =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex-1)         =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex)           = -4.0
                      poissonmatrix(localindex*columnsize+currentrowindex+1)         =  1.0
                  endif
              else
                  if(jj.EQ.1) then
                      poissonmatrix(localindex*columnsize+currentrowindex-numgrid_y) =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex)           = -4.0
                      poissonmatrix(localindex*columnsize+currentrowindex+1)         =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex+numgrid_y) =  1.0
                  else if(jj.EQ.numgrid_y) then
                      poissonmatrix(localindex*columnsize+currentrowindex-numgrid_y) =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex-1)         =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex)           = -4.0
                      poissonmatrix(localindex*columnsize+currentrowindex+numgrid_y) =  1.0
                  else
                      poissonmatrix(localindex*columnsize+currentrowindex-numgrid_y) =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex-1)         =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex)           = -4.0
                      poissonmatrix(localindex*columnsize+currentrowindex+1)         =  1.0
                      poissonmatrix(localindex*columnsize+currentrowindex+numgrid_y) =  1.0
                  endif
              endif
          enddo 
      enddo

      return

end subroutine construct_poissonmatrix
