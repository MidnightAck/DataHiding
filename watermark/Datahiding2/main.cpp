#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include"BMPop.h"
#include<fstream>
#include<vector>
#include <string.h>
#include "DistributionGen.h"

using namespace std;

#define a 100.8

vector<int> s;
vector<double> dct;

void showmat(vector<vector<double>> input) {
	for (int i = 0; i < input.size(); ++i) {
		for (int j = 0; j < input[0].size(); ++j) {
			cout << input[i][j] << "     ";
		}
		cout << endl;
	}
}

vector<vector<double>> DCT(vector<vector<double>> input, int row, int col)
{
	vector<vector<double>> output(row, vector<double>(col));
	//cout << "Test in DCT" << endl;
	double ALPHA, BETA;
	for (int u = 0; u < row; u++){
		for (int v = 0; v < col; v++){
			if (u == 0){
				ALPHA = sqrt(1.0 / row);
			}
			else{
				ALPHA = sqrt(2.0 / row);
			}

			if (v == 0){
				BETA = sqrt(1.0 / col);
			}
			else{
				BETA = sqrt(2.0 / col);
			}

			double tmp = 0.0;
			for (int i = 0; i < row; i++){
				for (int j = 0; j < col; j++){
					tmp += input[i][j] * cos((2 * i + 1)*u*PI / (2.0 * row)) * cos((2 * j + 1)*v*PI / (2.0 * col));
				}
			}
			output[u][v] = ALPHA * BETA * tmp;
		}
	}
	/*
	cout << "The result of DCT:" << endl;
	for (int i = 0; i < row; i++){
		for (int j = 0; j < col; j++){
			cout << output[i][j] << "     ";
		}
		cout << endl;
	}
	*/
	return output;
}

vector<vector<double>> IDCT(vector<vector<double>> input, int row, int col)
{
	vector<vector<double>> output(row, vector<double>(col));
	//cout << "Test in IDCT" << endl;
	double ALPHA, BETA;
	int u = 0;
	int v = 0;
	int i = 0;
	int j = 0;

	for (i = 0; i < row; i++)
	{
		for (j = 0; j < col; j++)
		{
			double tmp = 0.0;
			for (u = 0; u < row; u++)
			{
				for (v = 0; v < col; v++)
				{
					if (u == 0)
					{
						ALPHA = sqrt(1.0 / row);
					}
					else
					{
						ALPHA = sqrt(2.0 / row);
					}
					if (v == 0)
					{
						BETA = sqrt(1.0 / col);
					}
					else
					{
						BETA = sqrt(2.0 / col);
					}
					tmp += ALPHA * BETA * input[u][v] * cos((2 * i + 1)*u*PI / (2.0 * row)) * cos((2 * j + 1)*v*PI / (2.0 * col));
				}
			}
			
			output[i][j] = min(tmp,255);
			output[i][j] = max(tmp, 0);
			output[i][j] = floor(output[i][j]);
		}
	}

	return output;
}

void ASS(vector<int> w) {
	ifstream infile;
	infile.open("pixels.txt");

	ifstream infile1;
	infile1.open("logo.txt");
	vector<int> logo(1024);
	for (int i = 0; i < 1024; ++i)	infile1>>logo[i];

	ofstream outfile;
	outfile.open("DCT.txt");

	vector<vector<double>> bmp(10, vector<double>(10, 0));
	vector<vector<double>> middle(10, vector<double>(10, 0));
	for (int k = 0; k < 4096; ++k) {
		for (int i = 0; i < 8; ++i) {
			for (int j = 0; j < 8; ++j) {
				infile >> bmp[i][j];
			}
		}
		middle = DCT(bmp, 8, 8);
		//showmat(middle);
		double k1 = middle[3][4];
		double k2 = middle[4][3];
		
		double s1 = k1 + logo[k / 4] * (int)a * w[k * 2];
		double s2 = k2 + logo[k / 4] * (int)a * w[k * 2+1];
		s.push_back(s1);
		s.push_back(s2);

		dct.push_back(k1);
		dct.push_back(k2);

		middle[3][4] = s1;
		middle[4][3] = s2;

		bmp = IDCT(middle, 8, 8);
		for (int m = 0; m < 8; ++m) {
			for (int n = 0; n < 8; ++n) {
				outfile << bmp[m][n] << " ";
			}
			outfile << endl;
		}

	}
}

void decode(vector<int> w) {
	ifstream infile;
	infile.open("DCT.txt");

	ifstream infile1;
	infile1.open("logo.txt");
	vector<int> logo(1024);
	for (int i = 0; i < 1024; ++i)	infile1 >> logo[i];

	ofstream outfile;
	outfile.open("logo_decode.txt");

	vector<vector<double>> bmp(10, vector<double>(10, 0));
	vector<vector<double>> middle(10, vector<double>(10, 0));
	for (int k = 0; k < 4096; ++k) {
		for (int i = 0; i < 8; ++i) {
			for (int j = 0; j < 8; ++j) {
				infile >> bmp[i][j];
			}
		}
		middle = DCT(bmp, 8, 8);
		//showmat(middle);
		double k1 = middle[3][4];
		double k2 = middle[4][3];

		double s1 = k1 + logo[k / 4] * (int)a * w[k * 2];
		double s2 = k2 + logo[k / 4] * (int)a * w[k * 2 + 1];
		s.push_back(s1);
		s.push_back(s2);

		middle[3][4] = s1;
		middle[4][3] = s2;

		bmp = IDCT(middle, 8, 8);
		for (int m = 0; m < 8; ++m) {
			for (int n = 0; n < 8; ++n) {
				outfile << bmp[m][n] << " ";
			}
			outfile << endl;
		}

	}
}


void Randomize(vector<int>& w) {
	for (int i = 0; i < w.size(); ++i) {

		swap(w[i], w[Randint(i)]);
	}
}

int main()
{
	char BmpPath[] = "LENA.BMP";
	BMP2txt(BmpPath);
	char logoPath[] = "tj-logo.BMP";
	logo2txt(logoPath);

	vector<int> w(N, 1);
	for (int i = 500; i < N; ++i)	w[i] = -1;
	Randomize(w);

	ASS(w);
	
	txt2bmp();
	return 0;
}
