/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"

using namespace std;

double error(vector <double> v1, vector <double> v2, int n){
    if(n==0) return 0;
    else return sqrt((v1[n]-pow(v2[n],2))/n);
}

void data_blocking(vector <double> value, vector <double> &sum_prog, vector <double> &err_prog){
    vector <double> su2_prog;
    int N=sum_prog.size();
    
    for(int i=0;i<N;i++){
        sum_prog[i]=0;
        su2_prog.push_back(0);
    }
    
    for(int i=0;i<N;i++){
        for(int j=0;j<i+1;j++){
            sum_prog[i]=sum_prog[i]+value[j];
            su2_prog[i]=su2_prog[i]+pow(value[j],2);
        }
        sum_prog[i]=sum_prog[i]/(i+1);      //Cumulative average
        su2_prog[i]=su2_prog[i]/(i+1);      //Cumulative square average
        err_prog[i] = error(su2_prog,sum_prog,i);
    }
}

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
    
    //Parametri della simulazione
    double T=1.;
    double r=0.1;
    double sigma=0.25;
    double K=100;
    
    int N=100;
    int M=100000;
    int R=M/N;
    
    vector <double> sum_prog;
    vector <double> err_prog;
    
    vector <double> profit_call;
    vector <double> profit_put;
    
    for(int k=0;k<N;k++){
        sum_prog.push_back(0);
        err_prog.push_back(0);
        
        profit_call.push_back(0);
        profit_put.push_back(0);
    }
    
    ofstream output;
    output.open("output.txt");
    
    //DIRECT SAMPLING
    
    double S_0=100;
    double S=0;     //Value of S(T) with the direct sampling
    
    for(int i=0;i<N;i++){
        for(int j=0;j<R;j++){
            S=S_0*exp((r-0.5*pow(sigma,2))*T+sigma*rnd.Gauss(0.,T)*sqrt(T));
            profit_call[i]+=exp(-r*T)*max(0.,S-K);
            profit_put[i]+=exp(-r*T)*max(0.,K-S);
        }
        profit_call[i]/=R;     //call-option
        profit_put[i]/=R;      //put-option
    }
    
    //1) European call-option
    
    data_blocking(profit_call,sum_prog,err_prog);
    for(int i=0;i<N;i++){
        output<<sum_prog[i]<<" "<<err_prog[i]<<endl;
    }
    
    //2) European put-option
    
    data_blocking(profit_put,sum_prog,err_prog);
    for(int i=0;i<N;i++){
        output<<sum_prog[i]<<" "<<err_prog[i]<<endl;
    }
    
    //DISCRETIZED SAMPLING
    
    vector <double> S_d;    //Value of the asset price S(t_i) with i=0,...,99
    for(int i=0;i<N;i++){   //N=100=#(times in which [0,T] is divided)
        S_d.push_back(0);
        profit_call[i]=0;
        profit_put[i]=0;
    }
    S_d[0]=100;
    
    for(int i=1;i<N;i++){
        for(int j=0;j<R;j++){
            for(int k=1;k<N;k++){
                S_d[k]=S_d[k-1]*exp((r-0.5*pow(sigma,2))*(T/N)+sigma*rnd.Gauss(0,1)*sqrt(T/N));
            }
            profit_call[i]+=exp(-r*T)*max(0.,S_d[N-1]-K);
            profit_put[i]+=exp(-r*T)*max(0.,K-S_d[N-1]);
        }
        profit_call[i]/=R;     //call-option
        profit_put[i]/=R;      //put-option
    }
    
    //1) European call-option
    
    data_blocking(profit_call,sum_prog,err_prog);
    for(int i=0;i<N;i++){
        output<<sum_prog[i]<<" "<<err_prog[i]<<endl;
    }
    
    //2) European put-option
    
    data_blocking(profit_put,sum_prog,err_prog);
    for(int i=0;i<N;i++){
        output<<sum_prog[i]<<" "<<err_prog[i]<<endl;
    }
    
    output.close();
    
    rnd.SaveSeed();
    
    
    return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
