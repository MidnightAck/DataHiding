#include<iostream>
#include<fstream>

using namespace std;

#define N 10000



int main()
{
    for(int i=0;i<10;++i){
        cout<<GGDGen(1)<<endl;
    }
    cout<<endl<<endl;
    
    for(int i=0;i<10;++i){
        cout<<GGDGen(0.5)<<endl;
    }
    system("pause");
    return 0;
}
