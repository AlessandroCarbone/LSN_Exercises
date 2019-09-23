/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <vector>
#include <string>
#include "random.h"
#include "MolDyn_NVE.h"

using namespace std;

//Data blocking
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

int main(){ 
    Input();             //Inizialization
    int nconf = 1;
    int nblocks=200;
    for(int istep=1; istep <= nstep; ++istep){
        Move();           //Move particles with Verlet algorithm
        if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
        Measure();     //Properties measurement
        //ConfXYZ(nconf);   //Write actual configuration in XYZ format
        
        if(nconf==(nstep-1)){
            ConfFinal();
            system("mv config.final config.old");   //Write the second-last configuration to restart
        }
        nconf += 1;
    }
    
    ConfFinal();         //Write final configuration to restart
    Average(nblocks);
    
    return 0;
}


/**************************************************************************************************/

void Input(void){       //Prepare all stuff for the simulation
    ifstream ReadInput,ReadConf, ReadConfOld;

    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;

  
    //Random Numbers Generator
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
  
    ReadInput.open("input.dat");    //Read input
    
    ReadInput >> temp;
    
    ReadInput >> npart;
    cout << "Number of particles = " << npart << endl;
    
    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double)npart/rho;
    cout << "Volume of the simulation box = " << vol << endl;
    box = pow(vol,1.0/3.0);
    cout << "Edge of the simulation box = " << box << endl;
    
    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> iprint;
    ReadInput >> read;
    
    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of steps = " << nstep << endl << endl;
    ReadInput.close();
    
    //Prepare array for measurements
    iv = 0; //Potential energy
    ik = 1; //Kinetic energy
    ie = 2; //Total energy
    it = 3; //Temperature
    n_props = 4; //Number of observables
    
    double sumv2 = 0.0, fs;
    
    if(read==false){    //Just actual configuration
        
        //Qui vengono rinominati gli ultimi 2 files relativi alle configurazione finali della run precedente in modo da essere lette come files old.0 e old.final
        system("rm -f old.0");
        system("rm -f old.final");
        system("mv config.old old.final");
        system("mv config.final old.0");
        
        //Read initial configuration
        cout << "Read initial configuration from file config.0 " << endl << endl;
        ReadConf.open("config.0");
        for (int i=0; i<npart; ++i){
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] = Pbc(x[i] * box);
            y[i] = Pbc(y[i] * box);
            z[i] = Pbc(z[i] * box);
        }
        ReadConf.close();
    

        //Prepare initial velocities
        cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
        double sumv[3] = {0.0, 0.0, 0.0};
        for (int i=0; i<npart; ++i){
            vx[i] = rnd.Rannyu() - 0.5;
            vy[i] = rnd.Rannyu() - 0.5;
            vz[i] = rnd.Rannyu() - 0.5;

            sumv[0] += vx[i];
            sumv[1] += vy[i];
            sumv[2] += vz[i];
        }
    
        rnd.SaveSeed();
        for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
        
        for (int i=0; i<npart; ++i){
            vx[i] = vx[i] - sumv[0];
            vy[i] = vy[i] - sumv[1];
            vz[i] = vz[i] - sumv[2];

            sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        }
        sumv2 /= (double)npart;
        cout<<"Square velocity for particle = "<<sumv2<<endl<<endl;
        
        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
        for (int i=0; i<npart; ++i){
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;
            
            xold[i] = Pbc(x[i] - vx[i] * delta);
            yold[i] = Pbc(y[i] - vy[i] * delta);
            zold[i] = Pbc(z[i] - vz[i] * delta);
        }
    
    }
    
    if(read==true){         //Actual + Old configuration
        
        //The following vectors will be useful to contain the same information of the vectors xold, yold and zold because after the action of the Verlet algorithm those vectors will be overscripted with x,y and z.
        double a[npart], b[npart], c[npart];
        
        //Qui vengono rinominati gli ultimi 2 files relativi alle configurazione finali della run precedente in modo da essere lette come files old.0 e old.final
        system("rm -f old.0");
        system("rm -f old.final");
        system("mv config.old old.final");
        system("mv config.final old.0");
        
        //Read configuration old.0
        cout << "Read actual configuration from file old.0 " << endl << endl;
        ReadConf.open("old.0");
        for (int i=0; i<npart; ++i){
            ReadConf >> x[i] >> y[i] >> z[i];
            x[i] = Pbc(x[i] * box);
            y[i] = Pbc(y[i] * box);
            z[i] = Pbc(z[i] * box);
        }
        ReadConf.close();
        
        //Read configuration old.final
        cout << "Read old configuration from file old.final " << endl << endl;
        ReadConfOld.open("old.final");
        for (int i=0; i<npart; ++i){
            ReadConfOld >> xold[i] >> yold[i] >> zold[i];
            xold[i] = Pbc(xold[i] * box);
            yold[i] = Pbc(yold[i] * box);
            zold[i] = Pbc(zold[i] * box);
            
            a[i]=xold[i];
            b[i]=yold[i];
            c[i]=zold[i];
        }
        ReadConfOld.close();
        
        Move();     //One step of the Verlet algorithm
        
        for (int i=0; i<npart; ++i){
            vx[i]=Pbc(x[i]-a[i])/(2.0*delta);
            vy[i]=Pbc(y[i]-b[i])/(2.0*delta);
            vz[i]=Pbc(z[i]-c[i])/(2.0*delta);
            
            sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
        }
        
        sumv2 /= (double)(npart);
        cout<<"Square velocity for particle = "<<sumv2<<endl<<endl;
        
        fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
    
        for (int i=0; i<npart; ++i){
            vx[i] *= fs;
            vy[i] *= fs;
            vz[i] *= fs;
            
            a[i]=xold[i];
            b[i]=yold[i];
            c[i]=zold[i];
            
            xold[i] = Pbc(x[i] - 2.0 * vx[i] * delta);
            yold[i] = Pbc(y[i] - 2.0 * vy[i] * delta);
            zold[i] = Pbc(z[i] - 2.0 * vz[i] * delta);
            
            //I valori di xold, yold e zold trovati dall'algoritmo di Verlet corrispondono esattamente ai valori di x, y e z caricati dal file old.0
            x[i]=a[i];
            y[i]=b[i];
            z[i]=c[i];
        }
    }
    
    return;
}


void Move(void){    //Move particles with Verlet algorithm
    double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];
    
    for(int i=0; i<npart; ++i){     //Force acting on particle i
        fx[i] = Force(i,0);
        fy[i] = Force(i,1);
        fz[i] = Force(i,2);
    }

    for(int i=0; i<npart; ++i){     //Verlet integration scheme

        xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
        ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
        znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );
        
        vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
        vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
        vz[i] = Pbc(znew - zold[i])/(2.0 * delta);
        
        xold[i] = x[i];
        yold[i] = y[i];
        zold[i] = z[i];
        
        x[i] = xnew;
        y[i] = ynew;
        z[i] = znew;
    }
    return;
}

double Force(int ip, int idir){     //Compute forces as -Grad_ip V(r)
    double f=0.0;
    double dvec[3], dr;
    
    for (int i=0; i<npart; ++i){
        if(i != ip){
            dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
            dvec[1] = Pbc( y[ip] - y[i] );
            dvec[2] = Pbc( z[ip] - z[i] );
            
            dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
            dr = sqrt(dr);
            
            if(dr < rcut){
                f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8));       // -Grad_ip V(r)
            }
        }
    }
  
    return f;
}

void Measure(){     //Properties measurement
    //int bin;
    double v, t, vij, w, wij;
    double dx, dy, dz, dr;
    ofstream Epot, Ekin, Etot, Temp, Press;
    
    Epot.open("Results/output_epot.dat",ios::app);
    Ekin.open("Results/output_ekin.dat",ios::app);
    Temp.open("Results/output_temp.dat",ios::app);
    Etot.open("Results/output_etot.dat",ios::app);
    Press.open("Results/output_press.dat",ios::app);
    
    v = 0.0;    //reset observables
    t = 0.0;
    w = 0.0;
    
    //cycle over pairs of particles
    for (int i=0; i<npart-1; ++i){
        for (int j=i+1; j<npart; ++j){
            
            dx = Pbc( x[i] - x[j] );
            dy = Pbc( y[i] - y[j] );
            dz = Pbc( z[i] - z[j] );
            
            dr = dx*dx + dy*dy + dz*dz;
            dr = sqrt(dr);
            
            if(dr < rcut){
                //Potential energy
                vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
                v += vij;
                //Pressure
                wij=48.0/pow(dr,12) - 24.0/pow(dr,6);
                w+=wij;
            }
        }
    }
    
    //Kinetic energy
    for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    
    stima_pot = v/(double)npart; //Potential energy
    stima_kin = t/(double)npart; //Kinetic energy
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total enery
    stima_press=rho*stima_temp + w/(double)(3*vol);  //Pressure
    
    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Press << stima_press << endl;
    
    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Press.close();
    
    return;
}


void ConfFinal(void){   //Write final configuration
    ofstream WriteConf;
    
    cout << "Print final configuration" << endl << endl;
    WriteConf.open("config.final");
    
    for (int i=0; i<npart; ++i){
        WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    }
    WriteConf.close();
    return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
    ofstream WriteXYZ;
    
    WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
    WriteXYZ << npart << endl;
    WriteXYZ << "This is only a comment!" << endl;
    for (int i=0; i<npart; ++i){
        WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
    }
    WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

//Exercise_04.2
void Average(int nblocks){
    vector <double> ve,vu,vk,vt,vp,sum_prog,err_prog;
    
    for(int i=0;i<nblocks;i++){
        ve.push_back(0);
        vk.push_back(0);
        vu.push_back(0);
        vt.push_back(0);
        vp.push_back(0);
        
        sum_prog.push_back(0);
        err_prog.push_back(0);
    }
    
    ifstream input_e;
    ifstream input_k;
    ifstream input_u;
    ifstream input_t;
    ifstream input_p;
    
    input_e.open("Results/output_etot.dat");
    input_k.open("Results/output_ekin.dat");
    input_u.open("Results/output_epot.dat");
    input_t.open("Results/output_temp.dat");
    input_p.open("Results/output_press.dat");
    
    //Di seguito si memorizzano le varie misure di etot, ekin, epot e temp.
    double val=0;
    int L=nstep/nblocks;
    for(int i=0;i<nblocks;i++){
        for(int j=0;j<L;j++){
            input_e>>val;
            ve[i]+=val;
            input_k>>val;
            vk[i]+=val;
            input_u>>val;
            vu[i]+=val;
            input_t>>val;
            vt[i]+=val;
            input_p>>val;
            vp[i]+=val;
        }
        ve[i]/=L;
        vk[i]/=L;
        vu[i]/=L;
        vt[i]/=L;
        vp[i]/=L;
    }
    
    input_e.close();
    input_k.close();
    input_u.close();
    input_t.close();
    input_p.close();
    
    
    ofstream output_e;
    ofstream output_k;
    ofstream output_u;
    ofstream output_t;
    ofstream output_p;
    
    output_e.open("Results/ave_etot.dat");
    output_k.open("Results/ave_ekin.dat");
    output_u.open("Results/ave_epot.dat");
    output_t.open("Results/ave_temp.dat");
    output_p.open("Results/ave_press.dat");
    
    //Media dell'energia totale
    data_blocking(ve,sum_prog,err_prog);
    for(int k=0;k<nblocks;k++){
        output_e<<sum_prog[k]<<" "<<err_prog[k]<<endl;
    }
    
    //Media dell'energia cinetica
    data_blocking(vk,sum_prog,err_prog);
    for(int k=0;k<nblocks;k++){
        output_k<<sum_prog[k]<<" "<<err_prog[k]<<endl;
    }
    
    //Media dell'energia potenziale
    data_blocking(vu,sum_prog,err_prog);
    for(int k=0;k<nblocks;k++){
        output_u<<sum_prog[k]<<" "<<err_prog[k]<<endl;
    }
    
    //Media della temperatura
    data_blocking(vt,sum_prog,err_prog);
    for(int k=0;k<nblocks;k++){
        output_t<<sum_prog[k]<<" "<<err_prog[k]<<endl;
    }
    
    //Media della pressione (termine cinetico + termine potenziale)
    data_blocking(vp,sum_prog,err_prog);
    for(int k=0;k<nblocks;k++){
        output_p<<sum_prog[k]<<" "<<err_prog[k]<<endl;
    }
    
    output_e.close();
    output_k.close();
    output_u.close();
    output_t.close();
    output_p.close();
    
    return;
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
