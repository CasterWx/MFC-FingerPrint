#pragma once
#include<fstream>
#include<memory.h>
#include<string.h>
#include<iostream>
#include<atlimage.h>
#include<iomanip>

using namespace std;
// vs 编码为gbk
void Sort(unsigned char*data, int dsize) {
	unsigned char temp = 0;
	for (int i = 0; i < dsize; i++) {
		for (int j = dsize - 1; j > i; j--) {
			if (data[j] < data[j - 1]) {
				temp = data[j];
				data[j] = data[j - 1];
				data[j - 1] = temp;
			}
		}
	}
}


wchar_t*ToWideChar(char *str) {
	int num = MultiByteToWideChar(0, 0, str, -1, NULL, 0);
	wchar_t *wideStr = new wchar_t[num];
	MultiByteToWideChar(0, 0, str, -1, wideStr, num);
	return wideStr;
}

int ReadBMPImgFilePara(char* fileName, int &width, int &height, int &depth) {
	// 载入图像
	CImage image;
	HRESULT hResult = image.Load(ToWideChar(fileName));
	if (FAILED(hResult) || (image.IsNull())) {
		return -1;
	}
	// 获得图像参数
	width = image.GetWidth();
	height = image.GetHeight();
	depth = image.GetBPP();
	if (depth != 8) {
		return -2;
	}
	// 释放变量空间
	image.Destroy();

	return 0;
}

int ReadBMPImgFileData(char *fileName, unsigned char *data) {
	// 载入图像
	CImage image;
	HRESULT hResult = image.Load(ToWideChar(fileName));
	if (FAILED(hResult) || image.IsNull()) {
		return -1;
	}
	int width = image.GetWidth();
	int height = image.GetHeight();
	int depth = image.GetBPP();
	if (depth != 8) {
		return -2;
	}
	memset(data, 0, width*height);
	//读取图像数据
	int pitch = image.GetPitch();
	unsigned char *pData1 = (unsigned char*)image.GetBits();
	unsigned char* pData2 = data;
	unsigned char gray = 0;

	unsigned char *pRow1, *pRow2, *pPix1, *pPix2;
	for (int y = 0; y < height; y++) {
		pRow1 = pData1 + pitch * y;
		pRow2 = pData2 + width * y;
		for (int x = 0; x < width; x++) {
			pPix1 = pRow1 + x;
			gray = *pPix1;

			pPix2 = pRow2 + x;
			*pPix2 = gray;
		}
	}

	// 释放空间
	image.Destroy();
	return 0;
}
int WriteBMPImgFile(char *dstFileName, unsigned char** pusImgData) {
	FILE *fp = fopen(dstFileName, "r+b");
	if (!fp) {
		return -1;
	}

	int imgType, iWidth, iHeight;
	int iStartPos = 0;
	fseek(fp, 10L, SEEK_SET);
	fread((char*)(&iStartPos), 4, 1, fp);
	fseek(fp, 18L, SEEK_SET);
	fread((char*)(&iWidth), 4, 1, fp);
	fread((char*)(&iHeight), 4, 1, fp);
	unsigned short temp;
	fseek(fp, 28L, SEEK_SET);
	fread((char*)(&temp), 2, 1, fp);
	imgType = temp;
	if (imgType != 8) {
		return -2;
	}
	unsigned char* usImgData = *pusImgData;
	int iWidthInFile = 0;
	if (iWidth % 4 > 0) {
		iWidthInFile = iWidth - iWidth % 4 + 4;
	}
	else {
		iWidthInFile = iWidth;
	}
	for (int i = iHeight - 1; i >= 0; i--) {
		fseek(fp, iStartPos, SEEK_SET);
		fwrite((usImgData + i * iWidth), 1, iWidth, fp);
		iStartPos += iWidthInFile;
	}
	fclose(fp);
	return 0;
}
int SaveDataToImageFile(char* srcFile, char* dstFile, unsigned char* data) {
	CopyFile(ToWideChar(srcFile), ToWideChar(dstFile), false);
	WriteBMPImgFile(dstFile, &data);
	return 0;
}
int SaveDataToImageFile(char* srcFile, char* dstFile, float* data, float scale) {
	int iWidth, iHeight, iDepth;
	ReadBMPImgFilePara(srcFile, iWidth, iHeight, iDepth);
	CopyFile(ToWideChar(srcFile), ToWideChar(dstFile), false);
	unsigned char *tmpData = new unsigned char[iWidth*iHeight];
	for (int i = 0; i<int(iWidth*iHeight); i++) {
		tmpData[i] = unsigned char((scale* data[i]));
	}

	WriteBMPImgFile(dstFile, &tmpData);
	delete[] tmpData;
	return 0;
}
int SaveDataToTextFile(char* dstFile, unsigned char*data, int width, int height) {
	ofstream fout(dstFile, ios::out);
	if (!fout) {
		return -1;
	}
	// 按指定格式向文件写入数据
	int space = 5;
	for (int i = 0; i < height*width; i++) {
		fout << setw(space) << int(data[i]);
		if (i*width == (width - 1)) {
			fout << endl;
		}
	}
	fout.close();
	return 0;
}
int SaveDataToTextFile(char* dstFile, float *data, int width, int height) {
	ofstream fout(dstFile, ios::out);
	if (!fout) {
		return -1;
	}
	int preci = 6;
	int space = 16;
	fout.precision(preci);
	for (int i = 0; i < height*width - 1; i++) {
		fout << "  " << setw(space) << data[i];
		if (i%width == (width - 1)) {
			fout << endl;
		}
	}
	fout.close();
	return 0;
}
int ReadDatafromTextFile(char* srcFile, unsigned char* data, int iWidth, int iHeight) {
	// 打开源文件
	ifstream fin(srcFile, ios::in);
	if (!fin) {
		return -1;
	}
	// 读取数据
	int d = 0;
	for (int i = 0; i < iHeight*iWidth; i++) {
		fin >> d;
		data[i] = (unsigned char)d;
	}
	// 关闭文件
	fin.close();
	return 0;
}


int ReadDatafromTextFile(char* srcFile, float* data, int iWidth, int iHeight) {
	ifstream fin(srcFile, ios::in);
	if (!fin) {
		return -1;
	}
	for (int i = 0; i < iHeight*iWidth; i++) {
		fin >> data[i];
	}
	fin.close();
	return 0;
}


int ShowImageInCtrl(CStatic &picCtrl,char *filename) {
	CImage image;
	HRESULT hResult = image.Load(ToWideChar(filename));
	int width = image.GetWidth();	
	int height = image.GetHeight();

	CRect rect;
	picCtrl.GetClientRect(&rect);
	CDC *pDc = picCtrl.GetWindowDC();
	SetStretchBltMode(pDc->m_hDC,STRETCH_HALFTONE);
	image.StretchBlt(pDc->m_hAttribDC,rect,SRCCOPY);
	picCtrl.Invalidate(false);
	image.Destroy();
	picCtrl.ReleaseDC(pDc);
	return 0;
}
