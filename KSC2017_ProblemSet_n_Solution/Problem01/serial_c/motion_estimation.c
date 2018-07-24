#include "motion_estimation.h"

void DataProjectionH_C(unsigned char *src, int *nProjectionY, int nOffset, int nProjectionNum, int nWidth, int nHeight)
{
	int i, j;
	src += nOffset;

	memset(nProjectionY, 0, sizeof(int)*nHeight);

	for (j = 0; j < nHeight; j++)
	{
		for (i = 0; i < nProjectionNum; i++)
		{
			nProjectionY[j] += src[i + j*nWidth];
		}
	}
}

void DataProjectionV_C(unsigned char *src, int *nProjectionX, int nOffset, int nProjectionNum, int nWidth, int nHeight)
{
	int i, j;
	src += nOffset * nWidth;

	memset(nProjectionX, 0, sizeof(int)*nWidth);

	for (j = 0; j < nProjectionNum; j++)
	{
		for (i = 0; i < nWidth; i++)
		{
			nProjectionX[i] += src[i + j*nWidth];
		}
	}
}

int SAD_C(int *ref, int *cur, int nMV, int nOffset, int nMatchNum)
{
	int i, sad;

	ref += nOffset + nMV;
	cur += nOffset;

	sad = 0;
	for (i = 0; i<nMatchNum; i++)
	{
		sad += ABS(ref[i] - cur[i]);
	}

	return sad;
}

int LineMotionEstimation(int *nRef, int *nCurr, int nOffset, int nMatchNum, int nSearchRange)
{
	int i, min_sad=(1<<30), sad, nMV=0;

	for (i = -nSearchRange; i <= nSearchRange; i++)
	{
		sad = SAD_C(nRef, nCurr, i, nOffset, nMatchNum);
		if (sad < min_sad)
		{
			nMV = i;
			min_sad = sad;
		}
	}

	return nMV;
}

void MotionEstimation(unsigned char *pImage, int nWidth, int nHeight, int nSearchRange, int *nMV, int nFrameNum)
{
	int *nProjRef, *pProjCur;
	int i;
	int nMatchNumX, nMatchNumY;
	unsigned char *pRef, *pCur;

	nMatchNumX = nWidth - (nSearchRange << 1);
	nMatchNumY = nHeight - (nSearchRange << 1);

	nProjRef = (int*)malloc(sizeof(int)*nWidth);
	pProjCur = (int*)malloc(sizeof(int)*nWidth);

	for (i = 0; i < nFrameNum - 1; i++)
	{
		pRef = pImage + (long int)nWidth * (long int)nHeight * (long int)i;
		pCur = pRef + (long int)nWidth * (long int)nHeight;

		printf("Frame = %5d\n",i);

		DataProjectionV_C(pRef, nProjRef, nSearchRange, nMatchNumY, nWidth, nHeight);
		DataProjectionV_C(pCur, pProjCur, nSearchRange, nMatchNumY, nWidth, nHeight);
		nMV[2 * i + 0] = LineMotionEstimation(nProjRef, pProjCur, nSearchRange, nMatchNumX, nSearchRange);

		DataProjectionH_C(pRef, nProjRef, nSearchRange, nMatchNumX, nWidth, nHeight);
		DataProjectionH_C(pCur, pProjCur, nSearchRange + nMV[2 * i + 0], nMatchNumX, nWidth, nHeight);
		nMV[2 * i + 1] = LineMotionEstimation(nProjRef, pProjCur, nSearchRange, nMatchNumY, nSearchRange);

		DataProjectionV_C(pRef, nProjRef, nSearchRange, nMatchNumY, nWidth, nHeight);
		DataProjectionV_C(pCur, pProjCur, nSearchRange + nMV[2 * i + 1], nMatchNumY, nWidth, nHeight);
		nMV[2 * i + 0] = LineMotionEstimation(nProjRef, pProjCur, nSearchRange, nMatchNumX, nSearchRange);
	}

	free(nProjRef);
	free(pProjCur);
}
