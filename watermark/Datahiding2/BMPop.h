#define _CRT_SECURE_NO_WARNINGS
#include <math.h>
#include <iomanip> 
#include <stdlib.h>
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <malloc.h>

using namespace std;

unsigned char *pBmpBuf;//读入图像数据的指针

int bmpWidth;//图像的宽
int bmpHeight;//图像的高
RGBQUAD *pColorTable;//颜色表指针

int biBitCount;//图像类型，每像素位数


bool readBmp(char *bmpName){
	FILE *fp = fopen(bmpName, "rb");//二进制读方式打开指定的图像文件
	if (fp == 0)
		return 0;
	fseek(fp, sizeof(BITMAPFILEHEADER), 0);
	BITMAPINFOHEADER head;
	fread(&head, sizeof(BITMAPINFOHEADER), 1, fp); //获取图像宽、高、每像素所占位数等信息
	bmpWidth = head.biWidth;
	bmpHeight = head.biHeight;
	biBitCount = head.biBitCount;//定义变量，计算图像每行像素所占的字节数（必须是4的倍数）
	int lineByte = (bmpWidth * biBitCount / 8 + 3) / 4 * 4;//灰度图像有颜色表，且颜色表表项为256

	if (biBitCount == 8){
		//申请颜色表所需要的空间，读颜色表进内存
		pColorTable = new RGBQUAD[256];
		fread(pColorTable, sizeof(RGBQUAD), 256, fp);
	}

	pBmpBuf = new unsigned char[lineByte * bmpHeight];
	fread(pBmpBuf, 1, lineByte * bmpHeight, fp);
	fclose(fp);

	return 1;
}

//-----------------------------------------------------------------------------------------
//给定一个图像位图数据、宽、高、颜色表指针及每像素所占的位数等信息,将其写到指定文件中
bool saveBmp(char *bmpName, unsigned char *imgBuf, int width, int height, int biBitCount, RGBQUAD *pColorTable)
{
	if (!imgBuf)
		return 0;

	int colorTablesize = 0;
	if (biBitCount == 8)
		colorTablesize = 1024;

	int lineByte = (width * biBitCount / 8 + 3) / 4 * 4;
	FILE *fp = fopen(bmpName, "wb");
	if (fp == 0)
		return 0;
	//申请位图文件头结构变量，填写文件头信息
	BITMAPFILEHEADER fileHead;
	fileHead.bfType = 0x4D42;//bmp类型
	//bfSize是图像文件4个组成部分之和
	fileHead.bfSize = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER) + colorTablesize + lineByte * height;
	fileHead.bfReserved1 = 0;
	fileHead.bfReserved2 = 0;
	//bfOffBits是图像文件前3个部分所需空间之和
	fileHead.bfOffBits = 54 + colorTablesize;
	//写文件头进文件
	fwrite(&fileHead, sizeof(BITMAPFILEHEADER), 1, fp);
	//申请位图信息头结构变量，填写信息头信息
	BITMAPINFOHEADER head;
	head.biBitCount = biBitCount;
	head.biClrImportant = 0;
	head.biClrUsed = 0;
	head.biCompression = 0;
	head.biHeight = height;
	head.biPlanes = 1;
	head.biSize = 40;
	head.biSizeImage = lineByte * height;
	head.biWidth = width;
	head.biXPelsPerMeter = 0;
	head.biYPelsPerMeter = 0;

	fwrite(&head, sizeof(BITMAPINFOHEADER), 1, fp);
	if (biBitCount == 8)
		fwrite(pColorTable, sizeof(RGBQUAD), 256, fp);

	fwrite(imgBuf, height*lineByte, 1, fp);
	fclose(fp);
	return 1;
}

//----------------------------------------------------------------------------------------
//以下为像素的读取函数
void BMP2txt(char readPath[]) {
	readBmp(readPath);
	cout << "width=" << bmpWidth << " height=" << bmpHeight << " biBitCount=" << biBitCount << endl;

	int lineByte = (bmpWidth*biBitCount / 8 + 3) / 4 * 4;
	//cout << lineByte << endl;  512
	int m = 0, n = 0, count_xiang_su = 0;
	int count = 0;
	ofstream outfile("pixels.txt", ios::in | ios::trunc);
	for (int i = 0; i < bmpHeight; i += 8) {
		for (int j = 0; j < bmpWidth; j += 8) {
			for (int m = 0; m < 8; ++m) {
				for (int n = 0; n < 8; ++n) {
					outfile << (int)*(pBmpBuf + lineByte * (i + m) + j + n) << " ";
				}
				outfile << endl;
			}
		}
	}
	//将图像数据存盘

	char writePath[] = "nvcpy.BMP";
	saveBmp(writePath, pBmpBuf, bmpWidth, bmpHeight, biBitCount, pColorTable);
}

void txt2bmp() {
	ifstream infile;
	infile.open("DCT.txt");

	//ofstream out1;
	//out1.open("test.txt");
	unsigned char *pBmpBuf1;
	pBmpBuf1 = new unsigned char[512 * 512];
	int temp;
	int lineByte = (512 * biBitCount / 8 + 3) / 4 * 4;
	int m = 0, n = 0, count_xiang_su = 0;
	for (int i = 0; i < 512; i += 8) {
		for (int j = 0; j < 512; j += 8) {
			for (int m = 0; m < 8; ++m) {
				for (int n = 0; n < 8; ++n) {
					infile >> temp;
					*(pBmpBuf1 + lineByte * (i + m) + j + n) = temp;
					//out1 << (int)*(pBmpBuf1 + lineByte * (i + m) + j + n) << " ";
				}
				//out1 << endl;

			}
		}
	}

	//将图像数据存盘

	char writePath[] = "nvcpy1.BMP";
	saveBmp(writePath, pBmpBuf1, bmpWidth, bmpHeight, biBitCount, pColorTable);
}

unsigned char *read_1bmp(char *fname, int* _w, int* _h)
{
	unsigned char head[54];
	FILE *f = fopen(fname, "rb");

	// BMP header is 54 bytes
	fread(head, 1, 54, f);

	int w = head[18] + (((int)head[19]) << 8) + (((int)head[20]) << 16) + (((int)head[21]) << 24);
	int h = head[22] + (((int)head[23]) << 8) + (((int)head[24]) << 16) + (((int)head[25]) << 24);

	// lines are aligned on 4-byte boundary
	int lineSize = (w / 8 + (w / 8) % 4);
	int fileSize = lineSize * h;

	unsigned char *img = (unsigned char *)malloc(w * h), *data = (unsigned char *)malloc(fileSize);

	// skip the header
	fseek(f, 54, SEEK_SET);

	// skip palette - two rgb quads, 8 bytes
	fseek(f, 8, SEEK_CUR);

	// read data
	fread(data, 1, fileSize, f);

	// decode bits
	int i, j, k, rev_j;
	for (j = 0, rev_j = h - 1; j < h; j++, rev_j--) {
		for (i = 0; i < w / 8; i++) {
			int fpos = j * lineSize + i, pos = rev_j * w + i * 8;
			for (k = 0; k < 8; k++)
				img[pos + (7 - k)] = (data[fpos] >> k) & 1;
		}
	}

	free(data);
	*_w = w; *_h = h;
	return img;
}

void logo2txt(char readPath[]) {
	ofstream outfile;
	outfile.open("logo.txt");

	int w, h, i, j;
	unsigned char* img = read_1bmp(readPath, &w, &h);

	for (j = 0; j < h; j++)
	{
		for (i = 0; i < w; i++) {
			outfile << (img[j * w + i] ? 1 : -1) << " ";
			//cout << (img[j * w + i] ? '0' : '1') << " ";
		}
		//cout << endl;
		outfile << endl;
	}
}