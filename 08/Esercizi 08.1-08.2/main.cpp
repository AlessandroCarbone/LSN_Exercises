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

double trial_wave(double x, double mu, double sigma){
    return exp(-pow(x-mu,2)/(2*pow(sigma,2)))+exp(-pow(x+mu,2)/(2*pow(sigma,2)));
}

double potential(double x){
    return pow(x,4)-2.5*pow(x,2);
}

double D2_trial_wave(double x, double mu, double sigma){
    return exp(-pow(x+mu,2)/(2*pow(sigma,2)))*((1+exp(2*mu*x/pow(sigma,2)))*(pow(mu,2)-pow(sigma,2)+pow(x,2))-2*mu*x*(-1+exp(2*mu*x/pow(sigma,2))))/pow(sigma,4);
}

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
    
    double M=100000;
    double N=100;
    int L=M/N;
    
    //Exercise 08.1
    
    vector <double> H;
    for(int i=0;i<M;i++){
        H.push_back(0);
    }
    
    //Fixed parameters of the trial wave function
    double mu=1;
    double sigma=0.5;
    //The value of delta has been selected in order to get an acceptance rate of 50%
    double delta=0.95;
    double rap=0;
    
    double xold=1, xnew=0;
    H[0]=-0.5*D2_trial_wave(xold,mu,sigma)/trial_wave(xold,mu,sigma)+potential(xold);
    
    //Equilibration
    int N_eq=1000;
    for(int i=1;i<N_eq;i++){
        xnew=rnd.Gauss(xold,delta);
        H[i]=-0.5*D2_trial_wave(xnew,mu,sigma)/trial_wave(xnew,mu,sigma)+potential(xnew);
        rap=pow(trial_wave(xnew,mu,sigma),2)/pow(trial_wave(xold,mu,sigma),2);
        if(rap<1){
            if(rnd.Rannyu()>rap) H[i]=H[i-1];
            else xold=xnew;
        }
        else xold=xnew;
    }
    
    //Metropolis algorithm
    H[0]=-0.5*D2_trial_wave(xold,mu,sigma)/trial_wave(xold,mu,sigma)+potential(xold);
    int n_rej=0;
    for(int i=1;i<M;i++){
        xnew=rnd.Gauss(xold,delta);
        H[i]=-0.5*D2_trial_wave(xnew,mu,sigma)/trial_wave(xnew,mu,sigma)+potential(xnew);
        rap=pow(trial_wave(xnew,mu,sigma),2)/pow(trial_wave(xold,mu,sigma),2);
        if(rap<1){
            if(rnd.Rannyu()>rap){
                H[i]=H[i-1];
                n_rej++;
            }
            else xold=xnew;
        }
        else xold=xnew;
    }
    
    cout<<"Acceptance rate: "<<1-(double)n_rej/M<<endl;
    
    //Data blocking
    vector <double> ave;
    vector <double> sum_prog;
    vector <double> err_prog;
    
    for(int i=0;i<N;i++){
        ave.push_back(0);
        sum_prog.push_back(0);
        err_prog.push_back(0);
    }
    
    int k=0;
    for(int i=0;i<N;i++){
        for(int j=0;j<L;j++){
            k=j+i*L;
            ave[i]+=H[k];
        }
        ave[i]=ave[i]/L;
    }
    
    ofstream output_eigval;
    output_eigval.open("Results/ExpectedEigenvalue.txt");
    
    data_blocking(ave,sum_prog,err_prog);
    
    for(int i=0;i<N;i++){
        output_eigval<<sum_prog[i]<<" "<<err_prog[i]<<endl;
    }
    
    //Exercise 08.2
    
    double ground_state=0;
    double best_mu=1.;       //The most reliable mu associated with the ground state
    double best_sigma=0.5;  //The most reliable sigma associated with the ground state
    double sum=0;           //Analogue of sum_prog to be used only for the final average
    
    //Si utilizza l'ultimo xold del ciclo precedente
    
    //Cycle to evaluate the best parameters
    for(int i_par=0;i_par<400;i_par++){
        sigma=rnd.Gauss(best_sigma,0.5);
        mu=rnd.Gauss(best_mu,0.5);
        //cout<<xold<<endl;
        H[0]=-0.5*D2_trial_wave(xold,mu,sigma)/trial_wave(xold,mu,sigma)+potential(xold);
        
        for(int i=1;i<M;i++){
            xnew=rnd.Gauss(xold,delta);
            H[i]=-0.5*D2_trial_wave(xnew,mu,sigma)/trial_wave(xnew,mu,sigma)+potential(xnew);
            rap=pow(trial_wave(xnew,mu,sigma),2)/pow(trial_wave(xold,mu,sigma),2);
            if(rap<1){
                if(rnd.Rannyu()>rap) H[i]=H[i-1];
                else xold=xnew;
            }
            else xold=xnew;
        }
        
        sum=0;
        for(auto el: H){
            sum+=el;
        }
        sum/=M;
        
        if(sum<ground_state){
            ground_state=sum;
            //cout<<ground_state<<endl;
            best_sigma=sigma;
            best_mu=mu;
        }
    }
    
    cout<<"The approximated energy of the ground state is: E ="<<ground_state<<endl;
    cout<<"The parameters of the trial wave function are: mu = "<<best_mu<<" , sigma = "<<best_sigma<<endl;
    
    //Estimation of the ground state with {best_mu,best_sigma} and histogram for probability density
    
    xold=best_mu;
    H[0]=-0.5*D2_trial_wave(xold,best_mu,best_sigma)/trial_wave(xold,best_mu,best_sigma)+potential(xold);
    
    //Histogram for the pdf
    int nbins=100;
    double range=6.;        //pdf in range [-3.,3)
    vector <double> pdf;
    for(int i=0;i<nbins;i++){
        pdf.push_back(0.);
    }
    double bin_size=range/nbins;
    
    for(int i=1;i<M;i++){
        xnew=rnd.Gauss(xold,delta);
        H[i]=-0.5*D2_trial_wave(xnew,best_mu,best_sigma)/trial_wave(xnew,best_mu,best_sigma)+potential(xnew);
        rap=pow(trial_wave(xnew,best_mu,best_sigma),2)/pow(trial_wave(xold,best_mu,best_sigma),2);
        if(rap<1){
            if(rnd.Rannyu()>rap) H[i]=H[i-1];
            else xold=xnew;
        }
        else xold=xnew;
        //if(i<100) cout<<xold<<endl;
        
        //Filling the histogram
        for(int j=0;j<nbins;j++){
            if(xold>=(-3+j*bin_size) && xold<(-3+(j+1)*bin_size)){
                pdf[j]+=1.;
                break;
            }
        }
    }
    
    for(auto el: ave){
        el=0;
    }
    
    for(int i=0;i<N;i++){
        for(int j=0;j<L;j++){
            ave[i]+=H[j+i*L];
        }
        ave[i]/=L;
    }
    
    data_blocking(ave,sum_prog,err_prog);
    
    ofstream output_gs;
    output_gs.open("Results/GroundState.txt");

    for(int i=0;i<N;i++){
        output_gs<<sum_prog[i]<<" "<<err_prog[i]<<endl;
    }
    output_gs.close();
    
    ofstream output_pdf;
    output_pdf.open("Results/GroundState_pdf.txt");
    
    //We make the normalization of the pdf in the Python file
    for(int i=0;i<nbins;i++){
        output_pdf<<-3+i*bin_size<<" "<<pdf[i]<<endl;
    }
    output_pdf.close();
    
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
