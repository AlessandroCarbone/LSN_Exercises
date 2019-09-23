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
#include <vector>
#include <cmath>
//#define a_0 0.0529

using namespace std;

double error(vector <double> v1, vector <double> v2, int n){
    if(n==0) return 0;
    else return sqrt((v1[n]-pow(v2[n],2))/n);
}

double sph_arm1(double x,double y, double z){
    return exp(-2*sqrt(pow(x,2)+pow(y,2)+pow(z,2)))/M_PI;
}

double sph_arm2(double x,double y, double z){
    return pow(z,2)*exp(-sqrt(pow(x,2)+pow(y,2)+pow(z,2)))/(32*M_PI);
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
    int L=M/N;
    //double a_0=0.0529;
    
    vector <double> x;
    vector <double> y;
    vector <double> z;
    vector <double> r;      //Qui si inseriscono i valori di r
    
    vector <double> ave;
    vector <double> sum_prog;
    vector <double> err_prog;
    
    for(int i=0;i<M;i++){
        x.push_back(0);
        y.push_back(0);
        z.push_back(0);
        r.push_back(0);
    }
    
    for(int i=0;i<N;i++){
        ave.push_back(0);
        sum_prog.push_back(0);
        err_prog.push_back(0);
    }
    
    //Ground state wave function
    x[0]=0;
    y[0]=0;
    z[0]=0;
    r[0]=sqrt(pow(x[0],2)+pow(y[0],2)+pow(z[0],2));
    
    
    //Stima del delta migliore
    double delta=1;
    double rap=0;
    double n_rej=0;
    double val=0.9;     //Qui si inseriscono i valori di delta candidati
    double diff=1;
    for(int j=0;j<100;j++){
        x[0]=0;
        y[0]=0;
        z[0]=0;
        r[0]=sqrt(pow(x[0],2)+pow(y[0],2)+pow(z[0],2));
        n_rej=0;
        for(int i=1;i<10000;i++){
            x[i]=x[i-1]-val+2*val*rnd.Rannyu();
            y[i]=y[i-1]-val+2*val*rnd.Rannyu();
            z[i]=z[i-1]-val+2*val*rnd.Rannyu();
            r[i]=sqrt(pow(x[i],2)+pow(y[i],2)+pow(z[i],2));
            rap=sph_arm1(x[i],y[i],z[i])/sph_arm1(x[i-1],y[i-1],z[i-1]);
            if(rap<1){
                if(rnd.Rannyu()>rap){
                    x[i]=x[i-1];
                    y[i]=y[i-1];
                    z[i]=z[i-1];
                    r[i]=r[i-1];
                    n_rej++;
                }
            }
        }
        
        if(abs(1-n_rej/10000-0.5)<diff){
            delta=val;
            diff=abs(1-n_rej/10000-0.5);
        }
        val=val+0.005;
    }
    cout<<"Best choice for delta in sph_arm1"<<endl;
    cout<<delta<<" with acceptance probability equal to "<<0.5+diff<<endl;
    
    x[0]=0;
    y[0]=0;
    z[0]=0;
    r[0]=sqrt(pow(x[0],2)+pow(y[0],2)+pow(z[0],2));
    
    for(int i=1;i<M;i++){
        x[i]=x[i-1]-delta+2*delta*rnd.Rannyu();
        y[i]=y[i-1]-delta+2*delta*rnd.Rannyu();
        z[i]=z[i-1]-delta+2*delta*rnd.Rannyu();
        r[i]=sqrt(pow(x[i],2)+pow(y[i],2)+pow(z[i],2));
        rap=sph_arm1(x[i],y[i],z[i])/sph_arm1(x[i-1],y[i-1],z[i-1]);
        if(rap<1){
            if(rnd.Rannyu()>rap){
                x[i]=x[i-1];
                y[i]=y[i-1];
                z[i]=z[i-1];
                r[i]=r[i-1];
            }
        }
    }
    
    int k=0;
    for(int i=0;i<N;i++){
        for(int j=0;j<L;j++){
            k=j+i*L;
            ave[i]=ave[i]+r[k];
        }
        ave[i]/=L;
    }
    
    ofstream output1;
    output1.open("average.txt");
    ofstream output2;
    output2.open("points.txt");
    
    data_blocking(ave,sum_prog,err_prog);
    
    for(int k=0;k<M;k++){
        output2<<x[k]<<" "<<y[k]<<" "<<z[k]<<endl;
    }
    for(int k=0;k<N;k++){
        output1<<sum_prog[k]<<" "<<err_prog[k]<<endl;
    }
    
    //2p excited state wave function
    
    //Stima del delta migliore
    val=3;
    diff=1;
    for(int j=0;j<100;j++){
        x[0]=0;
        y[0]=0;
        z[0]=0;
        r[0]=sqrt(pow(x[0],2)+pow(y[0],2)+pow(z[0],2));
        n_rej=0;
        for(int i=1;i<10000;i++){
            x[i]=x[i-1]-val+2*val*rnd.Rannyu();
            y[i]=y[i-1]-val+2*val*rnd.Rannyu();
            z[i]=z[i-1]-val+2*val*rnd.Rannyu();
            r[i]=sqrt(pow(x[i],2)+pow(y[i],2)+pow(z[i],2));
            rap=sph_arm2(x[i],y[i],z[i])/sph_arm2(x[i-1],y[i-1],z[i-1]);
            if(rap<1){
                if(rnd.Rannyu()>rap){
                    x[i]=x[i-1];
                    y[i]=y[i-1];
                    z[i]=z[i-1];
                    r[i]=r[i-1];
                    n_rej++;
                }
            }
        }
        
        if(abs(1-n_rej/10000-0.5)<diff){
            delta=val;
            diff=abs(1-n_rej/10000-0.5);
        }
        val=val+0.005;
    }
    cout<<"Best choice for delta in sph_arm2"<<endl;
    cout<<delta<<" with acceptance probability equal to "<<0.5+diff<<endl;
    
    x[0]=0;
    y[0]=0;
    z[0]=1;
    r[0]=sqrt(pow(x[0],2)+pow(y[0],2)+pow(z[0],2));
    
    //delta=2;
    for(int i=1;i<M;i++){
        x[i]=x[i-1]-delta+2*delta*rnd.Rannyu();
        y[i]=y[i-1]-delta+2*delta*rnd.Rannyu();
        z[i]=z[i-1]-delta+2*delta*rnd.Rannyu();
        r[i]=sqrt(pow(x[i],2)+pow(y[i],2)+pow(z[i],2));
        rap=sph_arm2(x[i],y[i],z[i])/sph_arm2(x[i-1],y[i-1],z[i-1]);
        if(rap<1){
            if(rnd.Rannyu()>rap){
                x[i]=x[i-1];
                y[i]=y[i-1];
                z[i]=z[i-1];
                r[i]=r[i-1];
            }
        }
    }
    
    for(int i=0;i<N;i++){
        ave[i]=0;
    }
    
    for(int i=0;i<N;i++){
        for(int j=0;j<L;j++){
            k=j+i*L;
            ave[i]=ave[i]+r[k];
        }
        ave[i]/=L;
    }
    
    data_blocking(ave,sum_prog,err_prog);
    
    for(int k=0;k<M;k++){
        output2<<x[k]<<" "<<y[k]<<" "<<z[k]<<endl;
    }
    for(int k=0;k<N;k++){
        output1<<sum_prog[k]<<" "<<err_prog[k]<<endl;
    }
    
    output1.close();
    output2.close();
    
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
