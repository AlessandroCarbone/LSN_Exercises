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
#include <vector>
#include <cmath>
#include "random.h"

using namespace std;

double error(vector <double> v1, vector <double> v2, int n){
    if(n==0) return 0;
    else return pow((v1[n]-pow(v2[n],2))/n,0.5);
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
    
    int M= 1000000;
    int N= 100;
    int N_th=M/N;   //N_th=10000
    double d=5;     //Distance between the lines
    double L=1;     //Length of the needle
    
    //Si immagina che il bastoncino venga lanciato in uno spazio bidimensionale quadrato di area 100x100
    
    vector <double> ave;
    
    int N_hit=0;
    double x1,x2,y1,y2,mod,a,b;
    
    for(int i=0;i<N;i++){
        N_hit=0;
        for(int j=0;j<N_th;j++){
            x1=100*rnd.Rannyu();
            y1=100*rnd.Rannyu();
            mod=L+1;
            while(mod>L){
                a=(x1-L)+2*L*rnd.Rannyu();
                b=(y1-L)+2*L*rnd.Rannyu();
                mod=sqrt(pow(a-x1,2)+pow(b-y1,2));
            }
            
            if(b>=y1){
                x2=x1+L*(a-x1)/mod;
                y2=y1+L*sin(acos((a-x1)/mod));
            }
            else{
                x2=x1+L*cos(2*M_PI-acos((a-x1)/mod));
                y2=y1+L*sin(2*M_PI-acos((a-x1)/mod));
            }
            
            for(int k=0;k<=100;k=k+d){
                if(y1>k && y2<k) N_hit++;
                if(y2>k && y1<k) N_hit++;
                if(y1==k && y2==k) N_hit++;
            }
        }
        ave.push_back(2*L*N_th/(N_hit*d));  //Estimation of pi in each block
    }
    
    vector <double> sum_prog;
    vector <double> err_prog;
    
    for(int k=0;k<N;k++){
        sum_prog.push_back(0);
        err_prog.push_back(0);
    }
    
    data_blocking(ave,sum_prog,err_prog);
    
    ofstream output;
    output.open("output.txt");
    
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
