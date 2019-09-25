#include <iostream>
#include <fstream>
#include "Distribution.h"
using namespace std;

void menu(){
    system("cls");
    cout<<"================================"<<endl;
    cout<<"请选择功能："<<endl;
    cout<<"1.生成高斯分布"<<endl;
    cout<<"2.分析高斯分布参数"<<endl;
    cout<<"3.生成指数分布"<<endl;
    cout<<"4.分析指数分布参数"<<endl;
    cout<<"5.生成广义高斯分布"<<endl;
    cout<<"0.退出程序"<<endl;
    cout<<"================================"<<endl;
}

void GaussTest(){
    system("cls");
    DistributionGen DisGen;
    int n;
    cout<<"请输入数据规模： ";
    cin>>n;
    vector<double> dis=DisGen.GaussDisGen(n);
    ofstream outfile;
    outfile.open("GaussDis.txt");
    for(int 0;i<dis.size();i++) outfile<<dis[i]<<" ";
    outfile<<endl;
    cout<<"数据已输出到文件 GaussDis.txt"<<endl;
    while(1){
        cout<<"按1:将数据打印至控制台"<<endl;
        cout<<"按0:返回主菜单"<<endl;
        cin>>n;
        if(n==0)    return;
        if(n==1){
            for(int 0;i<dis.size();i++) cout<<dis[i]<<" ";
        }
        cout<<endl;
    }
    return;
}

void EXpTest(){
    system("cls");
    DistributionGen DisGen;
    int n;
    cout<<"请输入数据规模： ";
    cin>>n;
    double beta;
    cout<<"请输入指数（如果不输入默认为随机生成）"
    vector<double> dis=DisGen.ExpDisGen(n);
    ofstream outfile;
    outfile.open("ExpDis.txt");
    for(int 0;i<dis.size();i++) outfile<<dis[i]<<" ";
    outfile<<endl;
    cout<<"数据已输出到文件 ExpDis.txt"<<endl;
    while(1){
        cout<<"按1:将数据打印至控制台"<<endl;
        cout<<"按0:返回主菜单"<<endl;
        cin>>n;
        if(n==0)    return;
        if(n==1){
            for(int 0;i<dis.size();i++) cout<<dis[i]<<" ";
        }
        cout<<endl;
    }
    return;
}

int main()
{
    while(1){
        menu();
        int choice;
        cin>>choice;
        if(choice==0)   break;
        if(choice==1)   GaussTest();
        if(choice==2)   
    }
}
