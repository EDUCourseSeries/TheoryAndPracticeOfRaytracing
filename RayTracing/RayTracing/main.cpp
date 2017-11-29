#include <opencv2/opencv.hpp>
#include "Tracer.h"

#define IMAGE_WIDTH 800
#define IMAGE_HEIGHT 600

int main(void)
{

	CvSize csize;
	csize.width = IMAGE_WIDTH;
	csize.height = IMAGE_HEIGHT;

	IplImage *pImg = cvCreateImage(csize, 8, 3);

	CvScalar cs;
	for (int i = 0; i < IMAGE_HEIGHT; i++)
	{
		for (int j = 0; j < IMAGE_WIDTH; j++)
		{
			cs = cvGet2D(pImg, i, j);
			cs.val[2] = 255;
			cs.val[1] = 255;
			cs.val[0] = 255;
			cvSet2D(pImg, i, j, cs);
		}
	}

	cvShowImage("RayTracing", pImg);

	cvSaveImage("RayTracing.bmp", pImg);
	cvWaitKey();

	cvReleaseImage(&pImg);

	return 1;

}