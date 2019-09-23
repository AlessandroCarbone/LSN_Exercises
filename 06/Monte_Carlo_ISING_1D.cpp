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
#include <ostream>
#include <cmath>
#include <iomanip>
#include <string>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

//  Osservazioni
//  1) I files .0 nelle directories si riferiscono al caso metro=0 (Gibbs sampling), i files .1 al caso metro=1 (Metropolis sampling)
//  2) Le directory distinte con T1,...,T10 servono ad immagazzinare i dati riferiti alle varie temperature
//  3) input0.dat usa il Gibbs sampling e h=0.0, input1.dat usa il Metropolis sampling e h=0.0, input_mag1.dat il Metropolis sampling e h=0.02, input_mag0.dat il Metropolis sampling e h=0.02

int main(){
    double N=10;       //Numero di temperature con cui verificare l'andamento teorico
    double T=2;
    
    for(itemp=0;itemp<N;itemp++){
        
        //Creating directory where inserting data at fixed temperature
        string command="mkdir ";
        namedir="MCIsing_T"+to_string(itemp+1);   //Si contano le temperature prese con itemp
        command=command+namedir;
        system(command.c_str());
        
        h=0.0;
        while(h<0.03){
            metro=1;
            while(metro>-1){
                Input(T);   //Inizialization
                for(int iblk=1; iblk <= nblk; ++iblk){  //Simulation
                    Reset(iblk);   //Reset block averages
                    for(int istep=1; istep <= nstep; ++istep){
                        Move(metro);
                        Measure();
                        Accumulate();   //Update block averages
                    }
                    Averages(iblk);   //Print results for current block
                }
                ConfFinal(); //Write final configuration
                metro--;
            }
            h=h+0.02;
        }
        T=T-1.5/N;
    }
    return 0;
}


void Input(double T){
  
    ifstream ReadInput;

    cout << "Classic 1D Ising model             " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Nearest neighbour interaction      " << endl << endl;
    cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
    cout << "The program uses k_B=1 and mu_B=1 units " << endl;
    
    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();
    
    ifstream input("seed.in");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
    
    //Read input informations
    if(h==0.0){
        if(metro==0)ReadInput.open("input0.dat");
        else ReadInput.open("input1.dat");
    }
    if(h!=0.0){
        if(metro==0)ReadInput.open("input_mag0.dat");
        else ReadInput.open("input_mag1.dat");
    }
    
    ReadInput >> temp;
    if(T!=temp) temp=T;
    beta = 1.0/temp;
    cout << "Temperature = " << temp << endl;
    
    ReadInput >> nspin;
    cout << "Number of spins = " << nspin << endl;
    
    ReadInput >> J;
    cout << "Exchange interaction = " << J << endl;
    
    ReadInput >> h;
    cout << "External field = " << h << endl << endl;
    
    ReadInput >> metro; // if=1 Metropolis else Gibbs
    
    ReadInput >> nblk;
    
    ReadInput >> nstep;
    
    ReadInput >> read;
    
    if(metro==1) cout << "The program performs Metropolis moves" << endl;
    else cout << "The program performs Gibbs moves" << endl;
    
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;
    ReadInput.close();
    
    
    //Prepare arrays for measurements
    iu = 0;   //Energy
    iu2 = 1;  //Energy^2 to be used for the heat capacity
    is = 2;   //Sum of spins to be used for the susceptibility
    is2 = 3;  //(Sum of spins)^2 to be used for the magnetization
    
    n_props = 4; //Number of observables
    
    //Initial configuration
    
    //read è una variabile bool definita nel file .h che serve a stabilire se utilizzare la configurazione casuale o quella della precedente simulazione
    
    if(read==true){         //Random configuration
        for (int i=0; i<nspin; ++i){
            if(rnd.Rannyu() >= 0.5) s[i] = 1;
            else s[i] = -1;
        }
        system("rm -rf config.final");
    }
    else{               //Past configuration
        ifstream ReadConfig;
        ReadConfig.open("config.final");
        for(int i=0; i<nspin; ++i){
            ReadConfig >> s[i];
        }
        ReadConfig.close();
        system("rm -rf config.final");
        
    }
    
    
    //Evaluate energy etc. of the initial configuration
    Measure();
        
    //Print initial values for the potential energy and virial
    cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
    
}


void Move(int metro){
    int o;
    double energy_old, energy_new, sm;
    double alpha;
    double energy_up, energy_down;

    for(int i=0; i<nspin; ++i){
        //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
        o = (int)(rnd.Rannyu()*nspin);
        
        if(metro==1){     //Metropolis
            attempted++;
            sm=s[o];
            energy_old=Boltzmann(s[o],o);
            s[o]=-s[o];
            energy_new=Boltzmann(s[o],o);
            alpha=exp(-(energy_new-energy_old)*beta);
            if(alpha<1){
                if(rnd.Rannyu()>alpha) s[o]=sm;
            }
            if(s[o]==-sm) accepted++;
        }
        else{     //Gibbs sampling
            energy_up=Boltzmann(1,o);
            energy_down=Boltzmann(-1,o);
            alpha=1/(1+exp(-beta*(energy_down-energy_up)));
            if(rnd.Rannyu()<alpha) s[o]=1;
            else s[o]=-1;
        }
    }
}

double Boltzmann(int sm, int ip){
    double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
    return ene;
}

void Measure(){
    double u = 0.0, u2=0.0, sp=0.0, sp2=0.0;
    
    //cycle over spins
    for (int i=0; i<nspin; ++i){
        u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
        sp += s[i];
    }
    u2=pow(u,2);
    sp2=pow(sp,2);
    
    walker[iu] = u;
    walker[iu2] = u2;
    walker[is] = sp;
    walker[is2] = sp2;
}


void Reset(int iblk){ //Reset block averages

    if(iblk == 1){
       for(int i=0; i<n_props; ++i){
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i){
     blk_av[i] = 0;
   }
    
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void){  //Update block averages
    for(int i=0; i<n_props; ++i){
        blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
}


void Averages(int iblk){    //Print results for current block

    ofstream Ene, Heat, Mag, Chi;
    const int wd=12;
    
    cout << "Block number: " << iblk << endl;
    if(metro==1) cout << "Acceptance rate: " << accepted/attempted << endl << endl;
    
    
    if(h==0.0) Ene.open(namedir+"/output.ene."+to_string(metro),ios::app);
    else Ene.open(namedir+"/output.ene.mag."+to_string(metro),ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin;    //Energy
    glob_av[0]  += stima_u;
    glob_av2[0] += stima_u*stima_u;
    err_u=Error(glob_av[0],glob_av2[0],iblk);
    Ene << setw(wd) << iblk <<" "<<  setw(wd) << stima_u <<" "<< setw(wd) << glob_av[0]/(double)iblk << setw(wd) << err_u << endl;
    Ene.close();

    if(h==0.0)Heat.open(namedir+"/output.heatcapacity."+to_string(metro),ios::app);
    else Heat.open(namedir+"/output.heatcapacity.mag."+to_string(metro),ios::app);
    stima_c = pow(beta,2)*(blk_av[iu2]/blk_norm - pow(blk_av[iu]/blk_norm,2))/(double) nspin;   //Heat capacity
    glob_av[1]  += stima_c;
    glob_av2[1] += stima_c*stima_c;
    err_c=Error(glob_av[1],glob_av2[1],iblk);
    Heat << setw(wd) << iblk <<" "<<  setw(wd) << stima_c <<" "<< setw(wd) << glob_av[1]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();
    
    if(h==0.0) Chi.open(namedir+"/output.susceptibility."+to_string(metro),ios::app);
    else Chi.open(namedir+"/output.susceptibility.mag."+to_string(metro),ios::app);
    stima_x = beta*blk_av[is2]/blk_norm/(double) nspin;    //Susceptibility
    glob_av[2]  += stima_x;
    glob_av2[2] += stima_x*stima_x;
    err_x=Error(glob_av[2],glob_av2[2],iblk);
    Chi << setw(wd) << iblk <<" "<<  setw(wd) << stima_x <<" "<< setw(wd) << glob_av[2]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();
    
    if(h==0.0) Mag.open(namedir+"/output.magnetization."+to_string(metro),ios::app);
    else Mag.open(namedir+"/output.magnetization.mag."+to_string(metro),ios::app);
    stima_m = blk_av[is]/blk_norm/(double) nspin;     //Magnetization
    glob_av[3]  += stima_m;
    glob_av2[3] += stima_m*stima_m;
    err_m=Error(glob_av[3],glob_av2[3],iblk);
    Mag << setw(wd) << iblk <<" "<<  setw(wd) << stima_m <<" "<< setw(wd) << glob_av[3]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();
    
    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void){
    ofstream WriteConf;
    
    cout << "Print final configuration to file config.final " << endl << endl;
    WriteConf.open("config.final");
    for (int i=0; i<nspin; ++i){
        WriteConf << s[i] << endl;
    }
    WriteConf.close();
    
    rnd.SaveSeed();
}

int Pbc(int i){     //Algorithm for periodic boundary conditions
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk){
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
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
