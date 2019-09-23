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
        sum_prog[i]=sum_prog[i]/(i+1);
        su2_prog[i]=su2_prog[i]/(i+1);
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
    
   vector <double> r;
    int M= 100000;
    int N= 100;
    int L=M/N;
   for(int i=0;i<M;i++){
       r.push_back(rnd.Rannyu());
   }
    
// Punto 1
    vector <double> ave;
    vector <double> sum_prog;
    vector <double> err_prog;
    
    for(int l=0;l<N;l++){
        sum_prog.push_back(0);
        err_prog.push_back(0);
    }
    
    double sum=0;
    int k=0;
    
    for(int i=0;i<N;i++){
        sum=0;
        for(int j=0;j<L;j++){
            k=j+i*L;
            sum=sum+r[k];
        }
        ave.push_back(sum/L);
    }
    
    data_blocking(ave,sum_prog,err_prog);
    
    ofstream output;
    output.open("output.txt");
    
    for(int k=0;k<N;k++){
        output<<sum_prog[k]<<" "<<err_prog[k]<<endl;
    }
    
//Punto 2
    
    k=0;
    for(int i=0;i<N;i++){
        sum=0;
        for(int j=0;j<L;j++){
            k=j+i*L;
            sum=sum+pow(r[k]-0.5,2);
        }
        ave[i]=sum/L;
    }
    
    data_blocking(ave,sum_prog,err_prog);
    
    for(int k=0;k<N;k++){
        output<<sum_prog[k]<<" "<<err_prog[k]<<endl;
    }
    
// Parte 3
    
// Osservazione: Si sceglie di valutare i numeri contenuti nell'intervallo compreso tra 0 e 1/m.
    
    int m=100;
    vector <int> n;     //vector che conta i numeri in [0,1/m]
    vector <double> chi;
    for(int k=0;k<m;k++){
        n.push_back(0);
        chi.push_back(0);       //Questo perchè N=m=100
    }
    
    double rand=0;
    for(int i=0;i<N;i++){
        for(int l=0;l<m;l++){
            n[l]=0;
        }
        for(int j=0;j<10000;j++){
            rand=rnd.Rannyu();
            for(int k=0;k<m;k++){
                if(rand>=(double)k/m & rand<(double)(k+1)/m) n[k]++;
            }
        }
        for(int s=0;s<m;s++){
            chi[i]=chi[i]+pow(n[s]-10000/m,2)/(10000/m);
        }
    }
    
    data_blocking(chi,sum_prog,err_prog);
    
    for(int k=0;k<N;k++){
        output<<sum_prog[k]<<" "<<err_prog[k]<<endl;
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
