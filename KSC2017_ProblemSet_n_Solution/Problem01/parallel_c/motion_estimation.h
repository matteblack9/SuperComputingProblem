#ifndef __MOTION_ESTIMATION_H_
#define __MOTION_ESTIMATION_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define ABS(x)		(((x)>0) ? (x) : (-(x)))
#define MIN(x,y)	(((x)>(y)) ? (y) : (x))

void MotionEstimation(unsigned char *pImage, int nWidth, int nHeight, int nSearchRange, int *nMV, int nStaFrameNum, int nEndFrameNum, int myrank);

#endif

