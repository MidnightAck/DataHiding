#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "DistributionGen.h"
using namespace std;

#define GN 100000
#define a 1.8

void menu() {
	system("cls");
	cout << "================================" << endl;
	cout << "请选择功能：" << endl;
	cout << "1.生成高斯分布" << endl;
	cout << "2.分析高斯分布参数" << endl;
	cout << "3.生成指数分布" << endl;
	cout << "4.分析指数分布参数" << endl;
	cout << "5.生成广义高斯分布" << endl;
	cout << "0.退出程序" << endl;
	cout << "================================" << endl;
}

void GaussTest() {
	system("cls");
	DistributionGen DisGen;
	int n;
	cout << "请输入数据规模： ";
	cin >> n;
	vector<double> dis = DisGen.GaussDisGen(n);
	ofstream outfile;
	outfile.open("GaussDis.txt");
	for (int i = 0; i < n; i++) outfile << dis[i] << " ";
	outfile << endl;
	cout << "数据已输出到文件 GaussDis.txt" << endl;
	while (1) {
		cout << "按1:将数据打印至控制台" << endl;
		cout << "按0:返回主菜单" << endl;
		cin >> n;
		if (n == 0)    return;
		if (n == 1) {
			int count = 0;
			for (unsigned int i = 0; i < dis.size(); i++) {
				cout << setw(12) << dis[i];
				++count;
				if (count == 10) {
					count = 0;
					cout << endl;
				}
			}
		}
		cout << endl;
	}
	return;
}

void GaussAnalyzeTest() {
	system("cls");
	DistributionGen DisGen;
	vector<double> dis;
	ifstream infile;
	infile.open("GaussDis.txt");
	if (!infile) {
		cout << "请先生成高斯分布" << endl;
		while (1) {
			cout << "按0:返回主菜单" << endl;
			int n;
			cin >> n;
			if (n == 0)    return;
			cout << endl;
		}
		return;
	}
	for (int i = 0; i < 10000; i++) {
		double temp;
		infile >> temp;
		dis.push_back(temp);
	}
	cout << "已读入文件 GaussDis.txt" << endl;
	DisGen.AnalyzGaussDis(dis);
	int n;
	while (1) {
		cout << "按0:返回主菜单" << endl;
		cin >> n;
		if (n == 0)    return;
		cout << endl;
	}
	return;
}

void ExpTest() {
	system("cls");
	DistributionGen DisGen;
	int n;
	cout << "请输入数据规模： ";
	cin >> n;
	double beta;
	cout << "请输入指数（如果不输入默认为随机生成）";
	cin >> beta;
	vector<double> dis = DisGen.ExpDisGen(n, beta);
	ofstream outfile;
	outfile.open("ExpDis.txt");
	for (unsigned int i = 0; i < dis.size(); i++) outfile << dis[i] << " ";
	outfile << endl;
	cout << "数据已输出到文件 ExpDis.txt" << endl;
	while (1) {
		cout << "按1:将数据打印至控制台" << endl;
		cout << "按0:返回主菜单" << endl;
		cin >> n;
		if (n == 0)    return;
		if (n == 1) {
			int count = 0;
			for (unsigned int i = 0; i < dis.size(); i++) {
				cout << setw(12) << dis[i];
				++count;
				if (count == 10) {
					count = 0;
					cout << endl;
				}
			}
		}
		cout << endl;
	}
	return;
}

void ExpAnalyzeTest() {
	system("cls");
	DistributionGen DisGen;
	vector<double> dis;
	ifstream infile;
	infile.open("ExpDis.txt");
	if (!infile) {
		cout << "请先生成指数分布" << endl;
		while (1) {
			cout << "按0:返回主菜单" << endl;
			int n;
			cin >> n;
			if (n == 0)    return;
			cout << endl;
		}
		return;
	}
	for (int i = 0; i < 10000; i++) {
		double temp;
		infile >> temp;
		dis.push_back(temp);
	}
	cout << "已读入文件 ExpDis.txt" << endl;
	DisGen.AnalyzeExpDis(dis);
	int n;
	while (1) {
		cout << "按0:返回主菜单" << endl;
		cin >> n;
		if (n == 0)    return;
		cout << endl;
	}
	return;
}

vector<double> GGDTest(int n) {
	system("cls");
	DistributionGen DisGen;
	double beta = 1.0;
	vector<double> dis = DisGen.GGDDisGen(n, beta);
	//ofstream outfile;
	//outfile.open("GGDDis.txt");
	//for (unsigned int i = 0; i < dis.size(); i++) outfile << dis[i] << " ";
	//outfile << endl;
	//cout << "数据已输出到文件 GGDDis.txt" << endl;
	return dis;
}

void Randomize(vector<int>& w) {
	for (int i = 0; i < w.size(); ++i) {

		swap(w[i], w[Randint(i)]);
	}
}

double getPm() {

}



int main()
{
	vector<int> w(N,1);
	for (int i = 500; i < N; ++i)	w[i] = -1;
	Randomize(w);
	double ls = 0;
	
	for (int i = 0; i < GN; ++i) {
		vector<double> dis = GGDTest(1000);
		for (int j = 0; j < 1000; ++j) {
			ls += (dis[j] + a * w[j])*w[j] / N;
		}
		cout << i << endl;
	}
	ls /= GN;
	cout << ls << endl;
	return 0;
}
