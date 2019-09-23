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
#include "random.h"
#include <cmath>
#include <vector>

using namespace std;

// Parte 1
// Di seguito si vedono le distribuzioni di probabilità lorentziana ed esponenziale

double lorentz(double y, double mu, double gamma){
    return gamma*tan(M_PI*(y-0.5))+mu;
}

double expo(double y, double lambda){
    return -(1/lambda)*log(1-y);
}

double average(vector <double> v, int n){
    double sum=0;
    for(int i=0;i<n;i++){
        sum=sum+v[i];
    }
    return sum/n;
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


//Parte 2
    vector <double> s;  //vector che immagazzina i valori dello standard dice
    vector <double> e;  //vector che immagazzina i valori dell'exponential dice
    vector <double> l;  //vector che immagazzina i valori del lorentzian dice
    
    //Di seguito i vectors che immagazzinano il numero di volte che si ottiene come media un valore del dado tra 1 e 6.
    
    ofstream output;
    output.open("output.txt");
    
    // Caso N=1
    for(int i=0;i<10000;i++){
        s.push_back(rnd.Rannyu());
        e.push_back(expo(rnd.Rannyu(),1));
        l.push_back(lorentz(rnd.Rannyu(),0,1));
    }
    
    for(int k=0;k<10000;k++){
        output<<s[k]<<" "<<e[k]<<" "<<l[k]<<endl;
    }
    
    //Caso N=2

    //Nei vettori v1,v2,v3 si inseriscono le variabili random x_i di cui fare la media. Sapendo che bisogna fare la media anche per N=10 e N=100 si predispongono tali vettori in modo tale da avere dimensione uguale a 100.
    vector <double> v1,v2,v3;
    for(int i=0;i<100;i++){
        v1.push_back(0);
        v2.push_back(0);
        v3.push_back(0);
    }
    
    for(int i=0;i<10000;i++){
        for(int j=0;j<2;j++){
            v1[j]=rnd.Rannyu();
            v2[j]=expo(rnd.Rannyu(),1);
            v3[j]=lorentz(rnd.Rannyu(),0,1);
        }
        s[i]=average(v1,2);
        e[i]=average(v2,2);
        l[i]=average(v3,2);
    }
    
    for(int k=0;k<10000;k++){
        output<<s[k]<<" "<<e[k]<<" "<<l[k]<<endl;
    }
    
    //Caso N=10
    
    for(int i=0;i<10000;i++){
        for(int j=0;j<10;j++){
            v1[j]=rnd.Rannyu();
            v2[j]=expo(rnd.Rannyu(),1);
            v3[j]=lorentz(rnd.Rannyu(),0,1);
        }
        s[i]=average(v1,10);
        e[i]=average(v2,10);
        l[i]=average(v3,10);
    }
    
    for(int k=0;k<10000;k++){
        output<<s[k]<<" "<<e[k]<<" "<<l[k]<<endl;
    }
    
    //Caso N=100
    
    for(int i=0;i<10000;i++){
        for(int j=0;j<100;j++){
            v1[j]=rnd.Rannyu();
            v2[j]=expo(rnd.Rannyu(),1);
            v3[j]=lorentz(rnd.Rannyu(),0,1);
        }
        s[i]=average(v1,100);
        e[i]=average(v2,100);
        l[i]=average(v3,100);
    }
    
    for(int k=0;k<10000;k++){
        output<<s[k]<<" "<<e[k]<<" "<<l[k]<<endl;
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
