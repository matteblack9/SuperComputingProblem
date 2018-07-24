#!/usr/bin/env python24

#from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys



class PlotWave2D(object):
    def __init__(self, nx, ny, array_type):
        self.nx = nx
        self.ny = ny
        self.array_type = array_type

        #plt.ion()
        fig = plt.figure(figsize=(12,12))
        self.ax = fig.add_subplot(1,1,1)


        if array_type == 'F':
            self.plot_double_slit = self.plot_double_slit_y
        elif array_type == 'C':
            self.plot_double_slit = self.plot_double_slit_x


    def plot_double_slit_x(self, width, depth, gap):
        nx, ny = self.nx, self.ny
        x0 = nx/2
        ax = self.ax

        box1 = plt.Rectangle((x0, 0), depth, ny/2-(gap+width)/2, fc='k')
        box2 = plt.Rectangle((x0, ny/2-(gap-width)/2), depth, (gap-width), fc='k')
        box3 = plt.Rectangle((x0, ny/2+(gap+width)/2), depth, ny/2-(gap+width)/2, fc='k')
        ax.add_patch(box1)
        ax.add_patch(box2)
        ax.add_patch(box3)



    def plot_double_slit_y(self, width, depth, gap):
        nx, ny = self.nx, self.ny
        y0 = ny/2
        ax = self.ax

        box1 = plt.Rectangle((0, y0), nx/2-(gap+width)/2, depth, fc='k')
        box2 = plt.Rectangle((nx/2-(gap-width)/2, y0), (gap-width), depth, fc='k')
        box3 = plt.Rectangle((nx/2+(gap+width)/2, y0), nx/2-(gap+width)/2, depth, fc='k')
        ax.add_patch(box1)
        ax.add_patch(box2)
        ax.add_patch(box3)



    def read_binary_file(self, fpath):
        nx, ny = self.nx, self.ny
        array_type = self.array_type
		
        fp = open(fpath, 'rb')
        self.field = np.fromfile(fp, count=nx*ny, dtype=np.float64).reshape((nx,ny), order=array_type)



    def plot_field(self, vmin=None, vmax=None):
        ax = self.ax
        field = self.field

        img = ax.imshow(field.T, origin='lower', vmin=vmin, vmax=vmax)
        plt.colorbar(img)




if __name__ == '__main__':
    #-----------------------------------------
    # setup
    #-----------------------------------------
    '''
    #nx, ny = 2000, 2000
    nx, ny = 8000, 8000

    array_type = 'C'	# 'F':Fortran array, 'C': C array
    width, depth, gap = 120, 60, 1200
    '''
    #-----------------------------------------

    try:
        array_type = sys.argv[1]
        nx, ny = int(sys.argv[2]), int(sys.argv[3])

        if sys.argv[4] == 'slit':
            width, depth, gap = [int(x) for x in sys.argv[5:8]]
            fpath = sys.argv[8]
        else:
            fpath = sys.argv[4]

        pw = PlotWave2D(nx, ny, array_type)
        pw.read_binary_file(fpath)	
        pw.plot_field(vmin=-0.1, vmax=0.1)
        #pw.plot_field()

        if sys.argv[4] == 'slit':
            pw.plot_double_slit(width, depth, gap)

        plt.show()


    except IndexError:
        print 'Usage: ./plot_wave2d.py array_type nx ny [slit, width, depth, gap] field.bin'
        print 'example:'
        print '$ ./plot_wave2d.py F 2000 2000 field.bin'
        print '$ ./plot_wave2d.py F 8000 8000 field.bin'
        print '$ ./plot_wave2d.py F 2000 2000 slit 120 60 1200 field.bin'
        print '$ ./plot_wave2d.py F 8000 8000 slit 120 60 1200 field.bin'
        print ''
        print 'F: Fortran,  C: C'
        print ''
