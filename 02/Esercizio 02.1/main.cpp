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

    vector <double> r;
    int M= 100000;
    int N= 100;
    int L=M/N;
    for(int i=0;i<M;i++){
        r.push_back(rnd.Rannyu());
    }
    rnd.SaveSeed();
    
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
            sum=sum+(M_PI/2)*cos(M_PI*r[k]/2);
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
    //Per l'importance sampling si sceglie la funzione distribuita con probabilità p(x)=M_Pi/2*exp(-M_PI/2*x)/(1-exp(-M_PI/2)) (normalizzata nell'intervallo [0,1]).
    //Ho osservato che con questa scelta di p(x) il calcolo della varianza risulta essere più basso e ho mostrato il conto nel file di Mathematica contenuto nella cartella. Prima di selezionare questa funzione ho utilizzato più funzioni che allo stesso tempo avessero un andamento simile a quello della funzione da integrare e di cui fosse facile realizzare una distribuzione di probabilità. Un altro esempio di funzione p(x) che soddisfa questi 2 requisiti è la retta passante per i punti (0,M_PI/2) e (1,0).
    
    vector <double> p;
    for(int i=0;i<M;i++){
        p.push_back(-2/M_PI*log(1-(1-exp(-M_PI/2))*r[i]));
    }
    
    for(int i=0;i<N;i++){
        sum=0;
        for(int j=0;j<L;j++){
            k=j+i*L;
            sum=sum+(M_PI/2)*cos(M_PI*p[k]/2)/(M_PI/2*exp(-M_PI/2*p[k])/(1-exp(-M_PI/2)) );
        }
        ave[i]=(sum/L);
    }
    
    data_blocking(ave,sum_prog,err_prog);

    for(int k=0;k<N;k++){
        output<<sum_prog[k]<<" "<<err_prog[k]<<endl;
    }
    
    output.close();
    
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
