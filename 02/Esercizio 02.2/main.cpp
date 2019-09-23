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
    
    //Parte 1: CUBIC LATTICE
    int N_steps=100;
    int N=10000;
    
    vector <double> x;
    vector <double> y;
    vector <double> z;
    
    vector <double> r;
    vector <double> sum_prog;
    vector <double> err_prog;
    
    for(int k=0;k<N;k++){
        x.push_back(0);
        y.push_back(0);
        z.push_back(0);
    }
    
    for(int k=0;k<N_steps;k++){
        r.push_back(0);
        sum_prog.push_back(0);
        err_prog.push_back(0);
    }
    
    int sign=0;
    double rand=0;
    double val=0;
    
    for(int i=0;i<N;i++){
        for(int k=0;k<N;k++){
            x[k]=0;
            y[k]=0;
            z[k]=0;
        }
        val=0;
        for(int j=0;j<N_steps;j++){
            if(rnd.Rannyu()<0.5) sign=1;
            else sign=-1;
            rand=3*rnd.Rannyu();
            if(rand<1) x[j]=x[j]+sign;
            if(rand>1 && rand<2) y[j]=y[j]+sign;
            else z[j]=z[j]+sign;
            val=val+pow(x[j],2)+pow(y[j],2)+pow(z[j],2);
            r[j]=r[j]+val;
        }
    }
    
    for(int l=0;l<N_steps;l++){
        r[l]=sqrt(r[l]/N);
    }
    
    data_blocking(r,sum_prog,err_prog);
    
    ofstream output;
    output.open("output.txt");
    
    for(int i=0;i<N_steps;i++){
        output<<sum_prog[i]<<" "<<err_prog[i]<<endl;
    }
    
    
    //Parte 2: continuum
    for(int k=0;k<N_steps;k++){
        r[k]=0;
    }
    
    double theta=0;
    double phi=0;
    
    for(int i=0;i<N;i++){
        for(int k=0;k<N;k++){
            x[k]=0;
            y[k]=0;
            z[k]=0;
        }
        val=0;
        for(int j=0;j<N_steps;j++){
            theta=M_PI*rnd.Rannyu();
            phi=2*M_PI*rnd.Rannyu();
            
            x[j]=x[j]+cos(theta)*cos(phi);
            y[j]=y[j]+cos(theta)*sin(phi);
            z[j]=z[j]+sin(theta);
            val=val+pow(x[j],2)+pow(y[j],2)+pow(z[j],2);
            r[j]=r[j]+val;
        }
    }
    
    for(int l=0;l<N_steps;l++){
        r[l]=sqrt(r[l]/N);
    }
    
    data_blocking(r,sum_prog,err_prog);
    
    for(int i=0;i<N_steps;i++){
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
