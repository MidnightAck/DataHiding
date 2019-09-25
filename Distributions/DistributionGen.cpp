#include<iostream>
#include<fstream>

using namespace std;

#define N 10000

vector<double> gauss(N);
vector<double> cls(N);
//生成两个0-1之间的随机数，MT实现
double DistributionGen::RandGen(){
    std::random_device rd;  //获取随机数种子
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<double> dis(0.0, 1.0);
    
    double res=dis(gen);
    return res;
}

//生成两个0-1之间的随机数，MT实现
int DistributionGen::RandGenint(){
    std::random_device rd;  //获取随机数种子
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uniform_int_distribution<double> dis(0,10000);
    
    double res=dis(gen);
    return res;
}

//生成符合高斯分布的数据 Box-Muller
double DistributionGen::GaussGen(){
    //ofstream outfile;
    //outfile.open("Gauss.txt");
    vector<double> U(2);
    U[0]=RandGen();
    U[1]=RandGen();
    double Z=sqrt((-2)*log(U[0]))*sin(2*PI*U[1]);
    //outfile<<Z<<endl;
    //cout<<Z<<endl;gauss[i]=Z;
    
    //outfile.close();
    return Z;
}

//似然估计高斯分布的参数
void DistributionGen::AnalyzeGaussDis(){
    //ifstream infile;
    //infile.open("Gauss.txt");
    
    //vector<double> gauss(N);
    //for(int i=0;i<N;++i){
    //  infile>>gauss[i];
    //}
    double miu=0,sigma=0;
    for(int i=0;i<N;++i){
        miu+=gauss[i];
    }
    miu/=N;
    for(int i=0;i<N;++i){
        sigma+=(gauss[i]-miu)*(gauss[i]-miu);
    }
    sigma/=N;
    sigma=sqrt(sigma);
    cout<<miu<<endl<<sigma<<endl;
}

//生成指数分布z
double DistributionGen::ExpGen(double beita){
    // ofstream outfile;
    // outfile.open("clsRandom.txt");
    
    double pV=RandGen();
    double Z = (-1.0/beita)*log(1-pV);
    //outfile<<Z<<endl;
    cout<<Z<<endl;
    
    return Z;
}

void DistributionGen::analyzeExp(){
    
}

double DistributionGen::GammaGen(double alpha,double lambda)
{
    double u, v;
    double x, y;
    double b, c;
    double w, z;
    bool accept = false;
    
    if (alpha > 1.0){
        /* Best's rejection algorithm XG for gamma random variates (Best, 1978) */
        b = alpha - 1;
        c = 3*alpha-0.75;
        do {
            u = RandGen();
            v = RandGen();
            
            w = u*(1-u);
            y = sqrt(c/w)*(u-0.5);
            x = b+y;
            if (x>=0) {
                z = 64*w*w*w*v*v;
                double zz = 1-2*y*y/x;
                if (z-zz<1e-10){
                    accept = true;
                }
                else{
                    accept = false;
                }
                //accept=(z=1-2*y*y/x);
                //accept=(z==(1-2*y*y/x));
                if (!accept){
                    double logz = log(z);
                    double zzz = 2 * (b*log(x/b)-y);
                    if (logz-zzz<1e-10){
                        accept = true;
                    }
                    else{
                        accept = false;
                    }
                    //accept=(log(z)==2 * (b*log(x/b)-y));
                }
            }
        } while(!accept);
        return lambda*x;
    }
    else if (alpha == 1.0){
        return lambda * ExpGen(RandGen());
    }
    else if (alpha < 1.0){
        double dv = 0.0;
        double b=(exp(1.0)+alpha)/exp(1.0);
        do{
            double u1=RandGen();
            double p=b*RandGen();
            double y;
            if(p>1){
                y=-log((b-p)/alpha);
                if(u1<pow(y,alpha-1)){
                    dv=y;
                    break;
                }
            }else{
                y=pow(p,1/alpha);
                if(u1<exp(-y)){
                    dv=y;
                    break;
                }
            }
        }while(1);
        return dv*lambda;
    }
    return -1;
}

bool DistributionGen::BonGen(){
    double prob=0.5;
    default_random_engine generator;
    bernoulli_distribution b(prob);
    bool Z=0;
    int c=0;
    for(int i=0;i<RandGenint();++i){
        Z=b(generator);
        if(Z)   ++c;
    }
    cout<<c<<" ";
    return Z;
}

double GGDGen(double c){
    bool w=BonGen();
    double E=GammaGen(RandGen(),RandGen());
    double Z;
    if(w==1)    Z=pow(E,c);
    else    Z=-1*pow(E,c);
    return Z;
}

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
