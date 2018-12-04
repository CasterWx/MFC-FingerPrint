#include<fstream>
#include<memory.h>
#include<string.h>
#include<iostream>
#include<atlimage.h>
#include<iomanip>
#include"file.h"

using namespace std;

struct NEIGHBOR {
	int x;
	int y; 
	int type;
	float Theta;
	float Theta2Ridge;
	float Theta2ThisNibor;
	int distance;
};

struct MINUTIAE {
	int x;
	int y;
	int type;
	float theta;
	NEIGHBOR *neibors;
};

int ReadMinutiae(char* filename, MINUTIAE** minutiae) {
	FILE * fp = fopen(filename, "rb");
	if (!fp) {
		return -1;
	}
	const static int TemplateFileFlag = 0x3571027f;
	int flag;
	fread(&flag, sizeof(int), 1, fp);
	if (flag != TemplateFileFlag) {
		return -2;
	}
	int account;
	fread(&account, sizeof(int), 1, fp);
	*minutiae = new MINUTIAE[account];
	if (*minutiae == NULL) {
		return -3;
	}
	for (int i = 0; i < account; i++) {
		fread(&((*minutiae)[i]), sizeof(MINUTIAE), 1, fp);
	}

	fclose(fp);
	return account;
}

int MidFilter(unsigned char *ucImg, unsigned char *ucDstImg, int iWidth, int iHeight) {
	memset(ucDstImg, 0, iWidth*iHeight);
	unsigned char *pUp, *pDown, *pImg;
	unsigned char x[9];
	for (int i = 1; i < iHeight - 1; i++) {
		pUp = ucImg + (i - 1)*iWidth;
		pImg = ucImg + i * iWidth;
		pDown = ucImg + (i + 1)*iWidth;
		int j;
		for (j = 1; j < iWidth - 1; j++) {
			pUp++;
			pImg++;
			pDown++;
			x[0] = *(pUp - 1);
			x[1] = *(pImg - 1);
			x[2] = *(pDown - 1);
			x[3] = *pUp;
			x[4] = *pImg;
			x[5] = *pDown;
			x[6] = *(pUp + 1);
			x[7] = *(pImg + 1);
			x[8] = *(pDown + 1);

			Sort(x, 9);
			*(ucDstImg + i * iWidth + j) = x[4];
		}
	}
	int j;
	pDown = ucImg + iWidth;
	for (j = 1; j < iWidth - 1; j++) {
		x[0] = *(ucImg + j - 1);
		x[1] = *(ucImg + j);
		x[2] = *(ucImg + j + 1);
		x[3] = *(pDown + j - 1);
		x[4] = *(pDown + j);
		x[5] = *(pDown + j + 1);
		Sort(x, 6);
		*(ucDstImg + j) = x[3];
	}
	pUp = ucImg + iWidth * (iHeight - 2);
	pDown = ucImg + iWidth * (iHeight - 1);
	for (j = 1; j < iWidth - 1; j++) {
		x[0] = *(pDown + j - 1);
		x[1] = *(pDown + j);
		x[2] = *(pDown + j + 1);
		x[3] = *(pUp + j - 1);
		x[4] = *(pUp + j);
		x[5] = *(pUp + j + 1);
		Sort(x, 6);
		*(ucDstImg + iWidth * (iHeight - 1) + j) = x[3];
	}
	x[0] = *(ucImg);
	x[1] = *(ucImg + 1);
	x[2] = *(ucImg + iWidth);
	x[3] = *(ucImg + iWidth + 1);
	Sort(x, 4);
	*(ucDstImg) = x[2];

	x[0] = *(ucImg + iWidth - 1);
	x[1] = *(ucImg + iWidth - 2);
	x[2] = *(ucImg + 2 * iWidth - 1);
	x[3] = *(ucImg + 2 * iWidth - 2);

	Sort(x, 4);
	*(ucDstImg + iWidth - 1) = x[2];

	x[0] = *(ucImg + iWidth * (iHeight - 1));
	x[1] = *(ucImg + iWidth * (iHeight - 2));
	x[2] = *(ucImg + iWidth * (iHeight - 1) + 1);
	x[3] = *(ucImg + iWidth * (iHeight - 2) + 1);
	Sort(x, 4);
	*(ucDstImg + (iHeight - 1)*iWidth) = x[2];
	x[0] = *(ucImg + iWidth * (iHeight - 0) - 1);
	x[1] = *(ucImg + iWidth * (iHeight - 1) - 1);
	x[2] = *(ucImg + iWidth * (iHeight - 0) - 2);
	x[3] = *(ucImg + iWidth * (iHeight - 1) - 2);
	Sort(x, 4);

	*(ucDstImg + (iHeight - 0)*iWidth - 1) = x[2];
	 
	return 0;
}


int HistoNormalize(unsigned char* ucImg,unsigned char* ucNormImg, int iWidth,int iHeight) {
	unsigned int Histogram[256];

	memset(Histogram, 0, 256 * sizeof(int));
	for (int i = 0; i < iHeight; i++) {
		for (int j = 0; j < iWidth; j++) {
			Histogram[ucImg[i*iWidth + j]]++;
		}
	}
	double dMean = 0;
	for (int i = 1; i < 255; i++) {
		dMean += i * Histogram[i];
	}
	dMean = int(dMean/(iWidth*iHeight));
	double dSigma = 0;
	for (int i = 0; i < 255; i++) {
		dSigma += Histogram[i] * (i - dMean)*(i-dMean);
	}
	dSigma /= (iWidth*iHeight);
	dSigma = sqrt(dSigma);

	double dMean0 = 128, dSigma0 = 128;
	double dCoeff = dSigma0 / dSigma;
	for (int i = 0; i < iHeight; i++) {
		for(int j=0;j<iWidth;j++){
			double dVal = ucImg[i*iWidth + j]; 
			dVal = dMean0 + dCoeff * (dVal-dMean0);
			if (dVal < 0) {
				dVal = 0;
			}
			else if (dVal > 255) {
				dVal = 255;
			}
			ucNormImg[i*iWidth + j] = (unsigned char)dVal;
		}
	}
	return 0;
}


int ImgDirection(unsigned char* ucImg,float* fDirc,int iWidth,int iHeight) {
	const int SEMISIZ = 7;
	int dx[SEMISIZ * 2 + 1][SEMISIZ * 2 + 1];
	int dy[SEMISIZ * 2 + 1][SEMISIZ * 2 + 1];
	float fx, fy;
	memset(fDirc, 0, iWidth*iHeight * sizeof(float));
	for (int y = SEMISIZ + 1; y < iHeight - SEMISIZ - 1; y++) {
		for (int x = SEMISIZ + 1; x < iWidth - SEMISIZ - 1; x++) {
			for (int j = 0; j < SEMISIZ * 2 + 1; j++) {
				for (int i = 0; i < SEMISIZ * 2 + 1; i++) {
					int index1 = (y + j - SEMISIZ)*iWidth + x + i - SEMISIZ;
					int index2 = (y + j - SEMISIZ)*iWidth + x + i - SEMISIZ - 1;
					int index3 = (y + j - SEMISIZ - 1)*iWidth + x + i - SEMISIZ;
					dx[i][j] = int(ucImg[index1] - ucImg[index2]);
					dy[i][j] = int(ucImg[index1] - ucImg[index3]);
				}
			}
			fx = 0.0; 
			fy = 0.0;
			for (int j = 0; j < SEMISIZ * 2 + 1; j++) {
				for (int i = 0; i < SEMISIZ * 2 + 1; i++) {
					fx += 2 * dx[i][j] * dy[i][j];
					fy += (dx[i][j] * dx[i][j] - dy[i][j] * dy[i][j]);
				}
			}
			fDirc[y*iWidth + x] = atan2(fx,fy);
		}
	}
	return 0;
}


int DircLowPass(float *fDirc,float* fFitDirc,int iWidth,int iHeight) {
	const int DIR_FILTER_SIZE = 2;
	int blocksize = 2 * DIR_FILTER_SIZE + 1;
	int imgsize = iWidth * iHeight;
	float *filter = new float[blocksize*blocksize];
	float *phix = new float[imgsize];
	float *phiy = new float[imgsize];
	float *phi2x = new float[imgsize];
	float *phi2y = new float[imgsize];
	memset(fFitDirc,0,sizeof(float)*iWidth*iHeight);
	float tempSum = 0.0;
	for (int y = 0; y < blocksize; y++) {
		for (int x = 0; x < blocksize; x++) {
			filter[y*blocksize + x] = (float)(blocksize - (abs(DIR_FILTER_SIZE - x) + abs(DIR_FILTER_SIZE - y)));
			tempSum += filter[y*blocksize + x];
		}
	}
	for (int y = 0; y < blocksize; y++) {
		for (int x = 0; x < blocksize; x++) {
			filter[y*blocksize + x] /= tempSum;
		}
	}
	for(int y=0;y<iHeight;y++){
		for (int x = 0; x < iWidth; x++) {
			phix[y*iWidth + x] = cos(fDirc[y*iWidth + x]);
			phiy[y*iWidth + x] = sin(fDirc[y*iWidth + x]);
		}
	}
	memset(phi2x, 0, sizeof(float)*imgsize);
	memset(phi2y, 0, sizeof(float)*imgsize);
	float nx, ny;
	int val;
	for (int y = 0; y < iHeight - blocksize; y++) {
		for (int x = 0; x < iWidth - blocksize; x++) {
			nx = 0.0;
			ny = 0.0;
			for (int j = 0; j < blocksize; j++) {
				for (int i = 0; i < blocksize; i++) {
					val = (x + i) + (j + y)*iWidth;
					nx += filter[j*blocksize + i] * phix[val];
					ny += filter[j*blocksize + i] * phiy[val];
				}
			}
			val = x + y * iWidth;
			phi2x[val] = nx;
			phi2y[val] = ny;
		}
	}

	for (int y = 0; y < iHeight - blocksize; y++) {
		for (int x = 0; x < iWidth - blocksize; x++) {
			val = x + y * iWidth;
			fFitDirc[val] = atan2(phi2y[val],phi2x[val])*0.5;
		}
	}
	delete[] phi2y;
	delete[] phi2x;
	delete[] phiy;
	delete[] phix;

	return 0;
}


int Frequency(unsigned char* ucImg,float *fDirection,float* fFrequency,int iWidth,int iHeight){
	const int SIZE_L = 32;
	const int SIZE_W = 16;
	const int SIZE_L2 = 16;
	const int SIZE_W2 = 8;

	int peak_pos[SIZE_L];
	int peak_cnt;
	float peak_freq;
	float Xsig[SIZE_L];

	float dir = 0.0;
	float cosdir = 0.0;
	float sindir = 0.0;
	float maxPeak, minPeak;

	float *frequency1 = new float[iWidth*iHeight];
	memset(fFrequency, 0, sizeof(float)*iWidth*iHeight);
	memset(frequency1, 0, sizeof(float)*iWidth*iHeight);
	int x, y;
	int d, k;
	int u, v;
	for (y = SIZE_L2; y < iHeight - SIZE_L2; y++) {
		for (x = SIZE_L2; x < iWidth - SIZE_L2; x++) {
			dir = fDirection[(y + SIZE_W2)*iWidth + (x + SIZE_W2)];
			cosdir = -sin(dir);
			sindir = cos(dir);

			for (k = 0; k < SIZE_L; k++) {
				Xsig[k] = 0.0;
				for (d = 0; d < SIZE_W; d++) {
					u = (int)(x + (d - SIZE_W2)*cosdir + (k - SIZE_L2)*sindir);
					v = (int)(y + (d - SIZE_W2)*sindir - (k - SIZE_L2)*cosdir);
					if (u < 0) {
						u = 0;
					}
					else if (u > iWidth - 1) {
						u = iWidth - 1;
					}
					if (v < 0) {
						v = 0;
					}
					else if (v > iHeight - 1) {
						v = iHeight - 1;
					}
					Xsig[k] += ucImg[u + v * iWidth];
				}
				Xsig[k] /= SIZE_W;
			}
			maxPeak = minPeak = Xsig[0];
			for (k = 0; k < SIZE_L; k++) {
				if(minPeak>Xsig[k]){
					minPeak = Xsig[k];
				}
				if (maxPeak < Xsig[k]) {
					maxPeak = Xsig[k];
				}
			}
			peak_cnt = 0;
			if ((maxPeak - minPeak) > 64) {
				for(k=0;k<SIZE_L;k++){
					if ((Xsig[k - 1] < Xsig[k]) && (Xsig[k] >= Xsig[k + 1])) {
						peak_pos[peak_cnt++] = k;
					}
				}
			}
			peak_freq = 0.0;
			if (peak_cnt >= 2) {
				for (k = 0; k < peak_cnt - 1; k++) {
					peak_freq += (peak_pos[k + 1] - peak_pos[k]);
				}
				peak_freq /= peak_cnt - 1;
			}
			if (peak_freq<3.0 || peak_freq>25.0) {
				frequency1[x + y * iWidth] = 0.0;
			}
			else {
				frequency1[x + y * iWidth] = 1.0 / peak_freq;
			}
		}
	}
	for (y = SIZE_L2; y < iHeight - SIZE_L2; y++) {
		for (x = SIZE_L2; x < iWidth - SIZE_L2; x++) {
			k = x + y * iWidth;
			peak_freq = 0.0;
			for (v = -2; v <= 2; v++) {
				for (u = -2; u <= 2; u++) {
					peak_freq += frequency1[(x + u) + (y + v)*iWidth];
				}
			}
			fFrequency[k] = peak_freq/25;
		}
	}

	delete[] frequency1;
	return 0;
}


int GetMask(unsigned char* ucImg,float *fDirection,float *fFrequency,unsigned char *ucMask,int iWidth,int iHeight) {
	float freqMin = 1.0 / 25.0;
	float freqMax = 1.0 / 3.0;
	int x, y, k;
	int pos, posout;
	memset(ucMask,0,iWidth*iHeight);
	for (y = 0; y < iHeight; y++) {
		for (x = 0; x < iWidth; x++) {
			pos = x + y * iWidth;
			posout = x + y * iWidth;
			ucMask[posout] = 0;
			if (fFrequency[pos] >= freqMin && fFrequency[pos] <= freqMax) {
				ucMask[posout] = 255; 
			}
		}
	}
	for (k = 0; k < 4; k++) {
		for (y = 1; y < iHeight - 1; y++) {
			for (x = 1; x < iWidth - 1; x++) {
				if (ucMask[x + y * iWidth] == 0xFF) {
					ucMask[x - 1 + y * iWidth] |= 0x80;
					ucMask[x + 1 + y * iWidth] |= 0x80;
					ucMask[x + (y-1) * iWidth] |= 0x80;
					ucMask[x + (y+1) * iWidth] |= 0x80;
				}
			}
		}
		for (y = 1; y < iHeight - 1; y++) {
			for (x = 1; x < iWidth - 1; x++) {
				if (ucMask[x + y * iWidth]) {
					ucMask[x + y * iWidth] = 0xFF;
				}
			}
		}
	}
	for (k = 0; k < 12; k++) {
		for (y = 1; y < iHeight - 1; y++) {
			for (x = 1; x < iWidth - 1; x++) {
				if (ucMask[x + y * iWidth] == 0x0) {
					ucMask[x - 1 + y * iWidth] &= 0x80;
					ucMask[x + 1 + y * iWidth] &= 0x80;
					ucMask[x + (y - 1) * iWidth] &= 0x80;
					ucMask[x + (y + 1) * iWidth] &= 0x80;
				}
			}
		}
		for (y = 1; y < iHeight - 1; y++) {
			for (x = 1; x < iWidth - 1; x++) {
				if (ucMask[x + y * iWidth] != 0xFF) {
					ucMask[x + y * iWidth] = 0x0;
				}
			}
		}
	}
	return 0;
}


int GaborEnhance(unsigned char* ucImg, float* fDirection, float* fFrequency, unsigned char* ucMask, unsigned char* ucImgEnhanced, int iWidth, int iHeight) {
	const float PI = 3.141592654;
	int i, j, u, v;
	int wg2 = 5;
	float sum, f, g;
	float x2, y2;
	float dx2 = 1.0 / (4.0*4.0);
	float dy2 = 1.0 / (4.0*4.0);
	memset(ucImgEnhanced,0,iWidth*iHeight);
	for (j = wg2; j < iHeight - wg2; j++) {
		for (i = wg2; i < iWidth - wg2; i++) {
			if (ucMask[i + j * iWidth] == 0) {
				continue;
			}
			g = fDirection[i+j*iWidth];
			f = fFrequency[i+j*iWidth];
			g += PI / 2;
			sum = 0.0;
			for (v = -wg2; v <= wg2; v++) {
				for (u = -wg2; u <= wg2; u++) {
					x2 = -u * sin(g) + v * cos(g);
					y2 = u * cos(g) + v * sin(g);
					sum += exp(-0.5*(x2*x2*dx2 + y2 * y2*dy2))*cos(2 * PI*x2*f)*ucImg[(i - u) + (j - v)*iWidth];
				}
			}
			if (sum > 255.0) {
				sum = 255.0;
			}
			if (sum < 0.0) {
				sum = 0.0;
			}
			ucImgEnhanced[i + j * iWidth] = (unsigned char)sum;
		}
	}
	return 0;
}


int BinaryImg(unsigned char* ucImage,unsigned char* ucBinImage,int iWidth,int iHeight,unsigned char uThreshold) {
	unsigned char *pStart = ucImage, *pEnd = ucImage + iWidth * iHeight;
	unsigned char *pDest = ucBinImage;
	while (pStart < pEnd) {
		*pDest = *pStart > uThreshold ? 1 : 0;
		pStart++;
		pDest++;
	}
	return 0;
}


int BinaryToGray(unsigned char *ucBinImg,unsigned char *ucGrayImg,int iWidth,int iHeight) {
	unsigned char *pStart = ucBinImg, *pEnd = ucBinImg + iWidth * iHeight;
	unsigned char *pDest = ucGrayImg;

	while (pStart<pEnd) {
		*pDest = (*pStart) > 0 ? 255 : 0; 
		pStart++;
		pDest++;
	}
	return 0;
}


int Thinning(unsigned char *ucBinedImg,unsigned char *ucThinnedImage,int iWidth,int iHeight,int iIterativeLimit) {
	unsigned char x1, x2, x3, x4, x5, x6, x7, x8, xp;
	unsigned char g1, g2, g3, g4;
	unsigned char b1, b2, b3, b4;
	unsigned char np1, np2, npm;
	unsigned char *pUp, *pDown, *pImg;
	int iDeletePoints = 0;

	memcpy(ucThinnedImage,ucBinedImg,iWidth*iHeight);
	for (int it = 0; it < iIterativeLimit; it++) {
		iDeletePoints = 0;
		for (int i = 1; i < iHeight - 1; i++) {
			pUp = ucBinedImg + (i - 1)*iWidth;
			pImg = ucBinedImg + i * iWidth;
			pDown = ucBinedImg + (i + 1)*iWidth;
			for (int j = 1; j < iWidth - 1; j++) {
				pUp++;
				pImg++;
				pDown++;
				if (!*pImg) {
					continue;
				}
				x6 = *(pUp - 1);
				x5 = *(pImg - 1);
				x4 = *(pDown - 1);
				x7 = *pUp;
				xp = *pImg;
				x3 = *pDown;
				x8 = *(pUp + 1);
				x1 = *(pImg+1);
				x2 = *(pDown+1);


				b1 = !x1 && (x2 == 1 || x3 == 1);
				b2 = !x3 && (x4 == 1 || x5 == 1);
				b3 = !x5 && (x6 == 1 || x7 == 1);
				b4 = !x7 && (x8 == 1 || x1 == 1);

				g1 = (b1 + b2 + b3 + b4) == 1;

				np1 = x1 || x2;
				np1 += x3 || x4;
				np1 += x5 || x6;
				np1 += x7 || x8;
				np2  = x2 || x3;
				np2 += x4 || x5;
				np2 += x6 || x7;
				np2 += x8 || x1;

				npm = np1 > np2 ? np2 : np1;
				g2 = npm >= 2 && npm <= 3;
				g3 = (x1 && (x2 || x3 || !x8)) == 0;
				g4 = (x5 && (x6 || x7 || !x4)) == 0;

				if (g1&&g2&&g3) {
					ucThinnedImage[iWidth*i + j] = 0;
					++iDeletePoints;
				}
			}
		}
		memcpy(ucBinedImg,ucThinnedImage,iWidth*iHeight);
		for (int i = 1; i < iHeight - 1;i++) {
			pUp = ucBinedImg + (i - 1)*iWidth;
			pImg = ucBinedImg + i * iWidth;
			pDown = ucBinedImg + (i+1) * iWidth;
			for (int j = 1; j < iWidth - 1; j++) {
				pUp++;
				pImg++;
				pDown++;
				if (!*pImg) {
					continue;
				}
				x6 = *(pUp - 1);
				x5 = *(pImg - 1);
				x4 = *(pDown - 1);

				x7 = *pUp;
				xp = *pImg;
				x3 = *pDown;

				x8 = *(pUp + 1);
				x1 = *(pImg + 1);
				x2 = *(pDown + 1);

				b1 = !x1 && (x2 == 1 || x3 == 1);
				b2 = !x3 && (x4 == 1 || x5 == 1);
				b3 = !x5 && (x6 == 1 || x7 == 1);
				b4 = !x7 && (x8 == 1 || x1 == 1);

				g1 = (b1 + b2 + b3 + b4) == 1;

				np1 = x1 || x2;
				np1 += x3 || x4;
				np1 += x5 || x6;
				np1 += x7 || x8;

				np2 = x2 || x3;
				np2 += x4 || x5;
				np2 += x6 || x7;
				np2 += x8 || x1;

				npm = np1 > np2 ? np2 : np1;
				g2 = npm >= 2 && npm <= 3;

				g3 = (x1 && (x2 || x3 || !x8)) == 0;
				g4 = (x5 && (x6 || x7 || !x4)) == 0;

				if (g1&&g2&&g4) {
					ucThinnedImage[iWidth*i+j] = 0;
					++iDeletePoints;
				}
 			}
		}

		memcpy(ucBinedImg,ucThinnedImage,iWidth*iHeight);

		if (iDeletePoints == 0) {
			break;
		}
	}

	for (int i = 0; i < iHeight; i++) {
		for (int j = 0; j < iWidth; j++) {
			if (i < 16) {
				ucThinnedImage[i*iWidth + j] = 0;
			}
			else if (i >= iHeight - 16) {
				ucThinnedImage[i*iWidth + j] = 0;
			}
			else if (j < 16) {
				ucThinnedImage[i*iWidth + j] = 0;
			}
			else if (j >= iWidth - 16) {
				ucThinnedImage[i*iWidth + j] = 0;
			}
		}
	}
	return 0;
}


int Extract(unsigned char *ucThinImg,unsigned char *ucMinuImg,int iWidth,int iHeight){
	unsigned char *pDest = ucMinuImg;
	unsigned char *pUp, *pDown, *pImg;
	unsigned char x1, x2, x3, x4, x5, x6, x7, x8;
	unsigned char nc;
	int iMinuCount = 0;

	memset(pDest, 0, iWidth*iHeight);

	for (int i = 1; i < iHeight - 1; i++) {
		pUp = ucThinImg + (i - 1)*iWidth;
		pImg = ucThinImg + i * iWidth;
		pDown = ucThinImg + (i + 1)*iWidth;
		for (int j = 1; j < iWidth - 1; j++) {
			pUp++;
			pImg++;
			pDown++;
			if (!*pImg) {
				continue;
			}
			x6 = *(pUp - 1);
			x5 = *(pImg - 1);
			x4 = *(pDown - 1);
			x7 = *pUp;
			x3 = *pDown;
			x8 = *(pUp + 1);
			x1 = *(pImg + 1);
			x2 = *(pDown + 1);
			nc = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8;
			if (nc == 1) {
				pDest[i*iWidth + j] = 1;
				++iMinuCount;
			}
			else if (nc == 3) {
				pDest[i*iWidth + j] = 3;
				++iMinuCount;
			}
		}
	}
	return iMinuCount;
}


int CutEdge(MINUTIAE* minutiaes,int count,unsigned char*ucImg,int iWidth,int iHeight) {
	int minuCount = count;
	int x, y, type;
	bool del;
	int *pFlag = new int[minuCount];
	memset(pFlag,0,sizeof(int)*minuCount);
	for (int i = 0; i < minuCount; i++) {
		y = minutiaes[i].y - 1;
		x = minutiaes[i].x - 1;
		type = minutiaes[i].type ;
		del = true;
		if (x < iWidth / 2) {
			if (abs(iWidth / 2 - x) > abs(iHeight / 2 - y)) {
				while (--x >= 0) {
					if (ucImg[x + y * iWidth] > 0) {
						del = false;
						break;
					}
				}
			}
			else {
				if (y > iHeight / 2) {
					while (++y < iHeight) {
						if (ucImg[x + y * iWidth] > 0) {
							del = false;
							break;
						}
					}
				}
				else {
					while (--y == 0) {
						if (ucImg[x + y * iWidth] > 0) {
							del = false;
							break;
						}
					}
				}
			}
		}
		else {
			if (abs(iWidth / 2 - x) > abs(iHeight / 2 - y)) {
				while (++x < iWidth) {
					if (ucImg[x + y * iWidth] > 0) {
						del = false;
						break;
					}
				}
			}
			else {
				if (y > iHeight / 2) {
					while (++y < iHeight) {
						if (ucImg[x + y * iWidth] > 0) {
							del = false;
							break;
						}
					}
				}
				else {
					while (--y >= 0) {
						if (ucImg[x + y * iWidth] > 0) {
							del = false;
							break;
						}
					}
				}
			}
		}
		if (del) {
			pFlag[i] = 1;
			continue;
		}
	}
	int newCount = 0;
	for (int i = 0; i < minuCount; i++) {
		if (pFlag[i] == 0) {
			memcpy(&minutiaes[newCount],&minutiaes[i],sizeof(MINUTIAE));
			newCount++;
		}
	}
	delete[] pFlag;
	pFlag = NULL;
	return newCount;
}


int MinuFilter(unsigned char *minuData,unsigned char *thinData,MINUTIAE *minutiaes,int &minuCount,int iWidth,int iHeight) {
	float *dir = new float[iWidth*iHeight];
	memset(dir,0,iWidth*iHeight*sizeof(float));
	ImgDirection(thinData,dir,iWidth,iHeight);
	unsigned char* pImg;
	unsigned char val;
	int temp = 0;
	for (int i = 1; i < iHeight - 1; i++) {
		pImg = minuData + i * iWidth;
		for (int j = 1; j < iWidth - 1; j++) {
			++pImg;
			val = *pImg;
			if(val>0){
				minutiaes[temp].x = j + 1;
				minutiaes[temp].y = i + 1;
				minutiaes[temp].theta = dir[i*iWidth+j];
				minutiaes[temp].type = int(val);
				++temp;
			}
		}
	}
	delete[] dir; 
	minuCount = CutEdge(minutiaes,minuCount,thinData,iWidth,iHeight);
	int *pFlag = new int[minuCount];
	memset(pFlag,0,sizeof(int)*minuCount);
	int x1, x2, y1, y2, type1, type2;
	for (int i = 0; i < minuCount; i++) {
		x1 = minutiaes[i].x;
		y1 = minutiaes[i].y;
		type1 = minutiaes[i].type;
		for (int j = i + 1; j < minuCount; j++) {
			if (pFlag[i] == 1) {
				continue;
			}
			x2 = minutiaes[j].x;
			y2 = minutiaes[j].y;
			type2 = minutiaes[j].type;

			int r = (int)sqrt(float((y1-y2)*(y1-y2)+(x1-x2)*(x1-x2)));

			if (r <= 4) {
				if (type1 == type2) {
					if (type1 == 1) {
						pFlag[i] = pFlag[j] = 1;
					}
					else {
						pFlag[j] = 1;
					}
				}
				else if (type1 == 1) {
					pFlag[i] = 1;
				}
				else {
					pFlag[j] = 1;
				}
			}
		}

	}
	int newCount = 0;
	for (int i = 0; i < minuCount; i++) {
		if (pFlag[i] == 0) {
			memcpy(&minutiaes[newCount],&minutiaes[i],sizeof(MINUTIAE));
			newCount++;
		}
	}
	delete[] pFlag;
	minuCount = newCount;
	return 0;
}


int SaveMinutiae(MINUTIAE *minutiaes,int count,char *fileName) {
	FILE *fp = fopen(fileName,"wb");
	if (!fp) {
		return -1;
	}
	const static int TemplateFileFlag = 0x3571027f;
	fwrite(&TemplateFileFlag, sizeof(int), 1, fp);
	fwrite(&count,sizeof(int),1,fp);
	for (int i = 0; i < count; i++) {
		fwrite(&(minutiaes[i]),sizeof(MINUTIAE),1,fp);
	}
	fclose(fp);
	return 0;
}


float Angle2Points(int x1,int y1,int x2,int y2) {
	const float PI = 3.141592654;
	float diffY, diffX;
	float theta = 0.0;

	diffY = y2 - y1;
	diffX = x2 - x1;

	if (diffY < 0 && diffX>0) {
		theta = atan2(-1*diffY,diffX);
	}
	else if (diffY < 0 && diffX < 0) {
		theta = PI -  atan2(-1 * diffY, -1*diffX);
	}
	else if (diffY > 0 && diffX < 0) {
		theta = atan2(diffY, -1 * diffX);
	}
	else if (diffY > 0 && diffX > 0) {
		theta = PI - atan2(diffY, diffX);
	}
	else if (diffX == 0) {
		theta = PI / 2;
	}
	else {
		theta = 0.0;
	}
	return theta;
}


int BuildNabors(MINUTIAE *minutiae,int minuCount) {
	const int MAX_NEIGHBOR_EACH = 10;
	int x1, x2, y1, y2;
	int *pFlag = new int[minuCount];
	for (int i = 0; i < minuCount; i++) {
		x1 = minutiae[i].x;
		y1 = minutiae[i].y;

		memset(pFlag,0,sizeof(int)*minuCount);
		pFlag[i] = 1;
		minutiae[i].neibors = new NEIGHBOR[MAX_NEIGHBOR_EACH];
		if (minutiae[i].neibors == NULL) {
			return -1;
		}
		memset(minutiae[i].neibors,0,sizeof(NEIGHBOR)*MAX_NEIGHBOR_EACH);

		for (int neighborNo = 0; neighborNo < MAX_NEIGHBOR_EACH; neighborNo++) {
			int minDistance = 1000;
			int minNo = 0; 
			for (int j = 0; j < minuCount; j++) {
				if (pFlag[j] == 1) {
					continue;
				}
				x2 = minutiae[j].x;
				y2 = minutiae[j].y;
				int r = (int)sqrt(float((y1-y2)*(y1-y2)+(x1-x2)*(x1-x2)));

				if (r < minDistance) {
					minNo = j;
					minDistance = r;
				}
			}

			pFlag[minNo] = 1;
			minutiae[i].neibors[neighborNo].x = minutiae[minNo].x;
			minutiae[i].neibors[neighborNo].y = minutiae[minNo].y;
			minutiae[i].neibors[neighborNo].type = minutiae[minNo].type;
			minutiae[i].neibors[neighborNo].Theta = Angle2Points(minutiae[minNo].x,minutiae[minNo].y,x1,y1);

			minutiae[i].neibors[neighborNo].Theta2Ridge = minutiae[minNo].theta - minutiae[i].theta;

			minutiae[i].neibors[neighborNo].Theta2ThisNibor = minutiae[minNo].theta; 
			
			minutiae[i].neibors[neighborNo].distance = minDistance;

		}
	}
	delete[] pFlag;
	return 0;
}


float MinuSimilarity(int iWidth, int iHeight, MINUTIAE* minutiae1, int count1, MINUTIAE* minutiae2, int count2) {
	const int MAX_SIMILAR_PAIR = 100;
	const int MAX_NEIGHBOR_EACH = 10;
	BuildNabors(minutiae1, count1);
	BuildNabors(minutiae2, count2);

	int similarPair[MAX_SIMILAR_PAIR][2];

	memset(similarPair,0,100*2*sizeof(int));

	MINUTIAE *baseMinutiae;
	MINUTIAE *refMinutiae;

	int baseAccount, refAccount;

	if (count1 < count2) {
		baseMinutiae = minutiae1;
		baseAccount = count1;
		refMinutiae = minutiae2;
		refAccount = count2;
	}
	else {
		baseMinutiae = minutiae2;
		baseAccount = count2;
		refMinutiae = minutiae1;
		refAccount = count1;
	}

	NEIGHBOR *baseNeighbors = NULL;
	NEIGHBOR *refNeighbors = NULL;

	int similarMinutiae = 0;
	float baseTheta, refTheta;
	for (int i = 0; i < baseAccount; i++) {
		baseNeighbors = baseMinutiae[i].neibors;
		baseTheta = baseMinutiae[i].theta;
		int refSimilarNo = 0;
		int maxSimilarNerbors = 0;
		for (int j = 0; j < refAccount; j++) {
			if (refMinutiae[j].type != baseMinutiae[i].type) {
				continue;
			}
			refNeighbors = refMinutiae[j].neibors;
			refTheta = refMinutiae[j].theta;
			int thisSimilarNeigbors = 0;
			for (int m = 0; m < MAX_NEIGHBOR_EACH; m++) {
				for (int n = 0; n < MAX_NEIGHBOR_EACH; n++) {
					if (baseNeighbors[m].type != refNeighbors[n].type) {
						continue;
					}
					int dist = abs(int(baseNeighbors[m].distance-refNeighbors[n].distance));
					float theta1 = abs(float((baseNeighbors[m].Theta - baseTheta) - (refNeighbors[n].Theta - refTheta)));
					float theta2 = abs(float(baseNeighbors[m].Theta2Ridge - refNeighbors[n].Theta2Ridge));
					float theta3 = abs(float((baseNeighbors[m].Theta - baseNeighbors[m].Theta2ThisNibor) - (refNeighbors[n].Theta - refNeighbors[n].Theta2ThisNibor)));
					if (dist < 4 && theta1 < 0.15f&&theta2 < 0.15f&&theta3 < 0.15f) {
						++thisSimilarNeigbors;
						break;
					}
				}
			}
			if ((thisSimilarNeigbors >= MAX_NEIGHBOR_EACH * 3 / 10) && (similarMinutiae < MAX_SIMILAR_PAIR)) {
				similarPair[similarMinutiae][0] = i;
				similarPair[similarMinutiae][1] = refSimilarNo;
				++similarMinutiae;
			}
		}
	}

	float similarity = similarMinutiae/8.0f;
	similarity = similarMinutiae < 2 ? 0.0f : similarity;
	similarity = similarMinutiae > 8 ? 1.0f : similarity;
	return similarity;
}