/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "mpi.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "random.h"

using namespace std;

struct TSP_position{
    double x;
    double y;
};

double TSP_cost(vector <int> o,TSP_position p[]){
    double cost=0;
    int n=o.size();
    for(int i=1;i<n;i++){
        cost+=sqrt(pow(p[o[i]-1].x-p[o[i-1]-1].x,2)+pow(p[o[i]-1].y-p[o[i-1]-1].y,2));
    }
    cost+=sqrt(pow(p[o[0]-1].x-p[o[n-1]-1].x,2)+pow(p[o[0]-1].y-p[o[n-1]-1].y,2));
    return cost;
}

bool TSP_check(vector <int> o){
    int val=0;
    for(int i=0;i<o.size();i++){
        if(count(o.begin(),o.end(),i+1)>1){
            val=1;
            break;
        }
    }
    if(val==0) return true;     //Bonds respected
    else return false;          //Bonds not respected
}

double TSP_selection (vector <int> o1, vector <int> o2, TSP_position p[]){
    return TSP_cost(o2,p)/TSP_cost(o1,p);
}

//Mutations
void pair_permutation (vector <int> &vett, int a, int b){
    double box=vett[a];
    vett[a]=vett[b];
    vett[b]=box;
}

void shift (vector <int> &vett, int shift){
    vector <int> box;
    for(int i=0;i<vett.size();i++){
        box.push_back(vett[i]);
    }
    for(int i=0;i<vett.size()-shift;i++){
        vett[i+shift]=box[i];
    }
    for(int i=0;i<shift;i++){
        vett[i]=box[vett.size()-shift+i];
    }
}

//This operator works only if a+2*(b-a)<vett.size()
void shift (vector <int> &vett, int shift, int a, int b){   //a included, b excluded
    vector <int> box;
    for(int i=0;i<2*(b-a);i++){
        box.push_back(vett[a+i]);
    }
    ::shift (box,shift);
    
    for(int i=a;i<a+2*(b-a);i++){
        vett[i]=box[i];
    }
}

void invert_order(vector <int> &vett, int a, int b){     //a included,b included
    vector <int> box;
    for(int i=0;i<(b-a+1);i++){
        box.push_back(vett[b-i]);
    }
    for(int i=0;i<(b-a+1);i++){
        vett[a+i]=box[i];
    }
}

//This operator works only if a<b<c<d
void partial_permutation(vector <int> &vett, int a, int b, int c, int d){
    vector <int> box;
    
    for(int i=0;i<(b-a);i++){
        box.push_back(vett[a+i]);
    }
    for(int i=0;i<(d-c);i++){
        box.push_back(vett[c+i]);
    }
    
    random_shuffle(box.begin(),box.end());
    
    for(int i=0;i<(b-a);i++){
        vett[a+i]=box[i];
    }
    
    for(int i=0;i<(d-c);i++){
        vett[c+i]=box[b-a+i];
    }
}

void crossover(vector <int> &v1, vector <int> &v2, int cut){
    vector <int> box1, box2;
    for(int i=cut;i<v1.size();i++){
        box1.push_back(v1[i]);
        box2.push_back(v2[i]);
    }

    vector <int> pos;
    for(int i=0;i<box1.size();i++){
        for(int j=0;j<v2.size();j++){
            if(box1[i]==v2[j]){
                pos.push_back(j);
                break;
            }
        }
    }
    
    sort(pos.begin(),pos.end());
    for(int i=0;i<pos.size();i++){
        v1[cut+i]=v2[pos[i]];
    }
    
    for(auto el: pos) el=0;
    
    for(int i=0;i<box2.size();i++){
        for(int j=0;j<v1.size();j++){
            if(box2[i]==v1[j]){
                pos[i]=j;
                break;
            }
        }
    }
    
    sort(pos.begin(),pos.end());
    for(int i=0;i<pos.size();i++){
        v2[cut+i]=v1[pos[i]];
    }
}

double expo(double y, double lambda){       //To generate number exponentially distributed
    return -(1/lambda)*log(1-y);
}

//Osservazione: Per selezionare 2 vettori sulla base del valore di fitness minore (percorso più breve del venditore) si genera un numero con distribuzione esponenziale da sommare al valore minimo del fitness. Sucessivamente si sceglie il vettore che contiene il valore più vicino a quello ottenuto dalla somma.

//This function gives the index of the chromosome to use as a parent
int selection(vector <double> fitness, double rand){
    double min_fit=*min_element(fitness.begin(),fitness.end());
    double val=min_fit+expo(rand,15);
    int index=0;
    double min_diff=abs(fitness[0]-val);
    for(int k=1;k<fitness.size();k++){
        if(abs(fitness[k]-val)<min_diff){
            min_diff=abs(fitness[k]-val);
            index=k;
        }
    }
    return index;
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
        sum_prog[i]=sum_prog[i]/(i+1);
        su2_prog[i]=su2_prog[i]/(i+1);
        err_prog[i] = error(su2_prog,sum_prog,i);
    }
}

//The following struct to be used in MPI_Reduce
struct mpi_data{
    double l;
    int rank;
};

/////////////////////////////////////////////////////////
//MAIN
/////////////////////////////////////////////////////////

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
    
    int N=30;       //Number of cities
    
    //Here we create a vector which will be used to create the chromosomes
    vector <int> order;
    for(int i=0;i<N;i++){
        order.push_back(i+1);
    }
    
    TSP_position pos_c[N];
    
    //Punti su una circonferenza unitaria
    double angle=0;
    for(int i=0;i<N;i++){
        angle=2*M_PI*rnd.Rannyu();
        pos_c[i].x=cos(angle);
        pos_c[i].y=sin(angle);
    }
    
    ofstream circ_pos;
    circ_pos.open("Results/circ_position.txt");
    
    for(int i=0;i<N;i++){
        circ_pos<<pos_c[i].x<<" "<<pos_c[i].y<<endl;
    }
    
    TSP_position pos_s[N];
    
    double rand=0;
    
    //Punti dentro ad un quadrato
    for(int i=0;i<N;i++){
        rand=rnd.Rannyu();
        pos_s[i].x=rand;
        rand=rnd.Rannyu();
        pos_s[i].y=rand;
    }
    
    ofstream square_pos;
    square_pos.open("Results/square_position.txt");
    
    for(int i=0;i<N;i++){
        square_pos<<pos_s[i].x<<" "<<pos_s[i].y<<endl;
    }
    
    //Here we create two chromosomes which will be used for SA
    vector <int> chrom1 (N);
    vector <int> chrom2 (N);
    vector <int> new_chrom1 (N);
    vector <int> new_chrom2 (N);
    
    vector <int> v_rand;
    
    //MPI
    MPI::Init(argc,argv);
    int size = MPI::COMM_WORLD.Get_size();
    int rank = MPI::COMM_WORLD.Get_rank();
    
    for(int i=0;i<2*size;i++){
        v_rand.push_back(rnd.Rannyu(0,9999));
    }
    
    mpi_data irecv, isend;
    isend.rank=rank;
    
    /////////////////////////////////////////////////////////
    //CIRCUMFERENCE with R=1 and origin in (0,0)
    /////////////////////////////////////////////////////////
    
    random_shuffle(order.begin(),order.end());
    copy(order.begin(), order.end(), new_chrom1.begin());
    
    random_shuffle(order.begin(),order.end());
    copy(order.begin(), order.end(), new_chrom2.begin());
    
    rnd.SetRandom(seed,v_rand[2*rank],v_rand[2*rank+1]);
    
    double T=10;
    int steps=100;
    double min_cost=min(TSP_cost(new_chrom1,pos_c),TSP_cost(new_chrom2,pos_c));
    double cost_s=0, cost_p=0;
    
    vector <double> cost;
    
    for(int i=0;i<steps+1;i++){
        
        T=T-i*(double)T/steps;
        
        //Cycle of the i-th temperature
        for(int j=0;j<100;j++){
            
            copy(new_chrom1.begin(), new_chrom1.end(), chrom1.begin());
            copy(new_chrom2.begin(), new_chrom2.end(), chrom2.begin());
            
            //Crossover
            if(rnd.Rannyu()<0.7)
                crossover(new_chrom1, new_chrom2, (int)round(rnd.Rannyu(1,N)));
            
            //Mutation
            //Mutation of the second chromosome
            rand=rnd.Rannyu();
            if(rand<0.05) pair_permutation(new_chrom2,(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.05 && rand<0.1)
                shift(new_chrom2,(int)round(rnd.Rannyu(1,N)));
            if(rand>0.1 && rand<0.15) invert_order(new_chrom2,(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.15 && rand<0.2)
                partial_permutation(new_chrom2, (int)round(rnd.Rannyu(0,7)), (int)round(rnd.Rannyu(7,15)), (int)round(rnd.Rannyu(15,22)), (int)round(rnd.Rannyu(22,30)));
            
            //Mutation of the first chromosome
            rand=rnd.Rannyu();
            if(rand<0.05) pair_permutation(new_chrom1,(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.05 && rand<0.1)
                shift(new_chrom1,(int)round(rnd.Rannyu(1,N)));
            if(rand>0.1 && rand<0.15) invert_order(new_chrom1,(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.15 && rand<0.2)
                partial_permutation(new_chrom1, (int)round(rnd.Rannyu(0,7)), (int)round(rnd.Rannyu(7,15)), (int)round(rnd.Rannyu(15,22)), (int)round(rnd.Rannyu(22,30)));
            
            //Metropolis algorithm to obtain SA
            cost_p=min(TSP_cost(chrom1,pos_c),TSP_cost(chrom2,pos_c));
            
            //Second chromosome
            cost_s=TSP_cost(new_chrom2,pos_c);
            
            if(cost_s>cost_p){
                if(rnd.Rannyu()>exp(-(cost_s-cost_p)/T))
                    copy(chrom2.begin(), chrom2.end(), new_chrom2.begin());
            }
            
            //First chromosome
            cost_s=TSP_cost(new_chrom1,pos_c);
            if(cost_s>cost_p){
                if(rnd.Rannyu()>exp(-(cost_s-cost_p)/T))
                    copy(chrom1.begin(), chrom1.end(), new_chrom1.begin());
            }
            
            cost_s=min(TSP_cost(new_chrom1,pos_c),TSP_cost(new_chrom2,pos_c));
            
            if(cost_s<min_cost)
                min_cost=min(TSP_cost(new_chrom1,pos_c),TSP_cost(new_chrom2,pos_c));
            
        }
        
        cost.push_back(min_cost);
        
    }
    
    isend.l=min(TSP_cost(new_chrom1,pos_c),TSP_cost(new_chrom2,pos_c));
    cout<<"Circumference - Rank "<<isend.rank<<" - Cost : "<<isend.l<<endl;
    MPI_Reduce(&isend, &irecv, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI::COMM_WORLD);
    if(rank==0)
        cout<<"Rank "<<irecv.rank<<" has the lowest value: Cost = "<<irecv.l<<endl;
    
    circ_pos.close();
    
    /////////////////////////////////////////////////////////
    //SQUARE with L=1: x in (0,1), y in (0,1)
    /////////////////////////////////////////////////////////
    
    //We initialize the next cycles using the chromosomes used at the end of the circumference algorithm
    T=4;
    for(int i=0;i<steps;i++){
        
        T=T-i*(double)T/steps;
        
        //Cycle of the i-th temperature
        for(int j=0;j<100;j++){
            
            copy(new_chrom1.begin(), new_chrom1.end(), chrom1.begin());
            copy(new_chrom2.begin(), new_chrom2.end(), chrom2.begin());
            
            //Crossover
            if(rnd.Rannyu()<0.7)
                crossover(new_chrom1, new_chrom2, (int)round(rnd.Rannyu(1,N)));
            
            //Mutation
            //Mutation of the second chromosome
            rand=rnd.Rannyu();
            if(rand<0.05) pair_permutation(new_chrom2,(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.05 && rand<0.1)
                shift(new_chrom2,(int)round(rnd.Rannyu(1,N)));
            if(rand>0.1 && rand<0.15) invert_order(new_chrom2,(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.15 && rand<0.2)
                partial_permutation(new_chrom2, (int)round(rnd.Rannyu(0,7)), (int)round(rnd.Rannyu(7,15)), (int)round(rnd.Rannyu(15,22)), (int)round(rnd.Rannyu(22,30)));
            
            //Mutation of the first chromosome
            rand=rnd.Rannyu();
            if(rand<0.05) pair_permutation(new_chrom1,(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.05 && rand<0.1)
                shift(new_chrom1,(int)round(rnd.Rannyu(1,N)));
            if(rand>0.1 && rand<0.15) invert_order(new_chrom1,(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.15 && rand<0.2)
                partial_permutation(new_chrom1, (int)round(rnd.Rannyu(0,7)), (int)round(rnd.Rannyu(7,15)), (int)round(rnd.Rannyu(15,22)), (int)round(rnd.Rannyu(22,30)));
            
            //Metropolis algorithm to obtain SA
            cost_p=min(TSP_cost(chrom1,pos_s),TSP_cost(chrom2,pos_s));
            
            //Second chromosome
            cost_s=TSP_cost(new_chrom2,pos_s);
            
            if(cost_s>cost_p){
                if(rnd.Rannyu()>exp(-(cost_s-cost_p)/T))
                    copy(chrom2.begin(), chrom2.end(), new_chrom2.begin());
            }
            
            //First chromosome
            cost_s=TSP_cost(new_chrom1,pos_s);
            if(cost_s>cost_p){
                if(rnd.Rannyu()>exp(-(cost_s-cost_p)/T))
                    copy(chrom1.begin(), chrom1.end(), new_chrom1.begin());
            }
            
            cost_s=min(TSP_cost(new_chrom1,pos_s),TSP_cost(new_chrom2,pos_s));
            
            if(cost_s<min_cost)
                min_cost=min(TSP_cost(new_chrom1,pos_s),TSP_cost(new_chrom2,pos_s));
            
        }
        
        cost[i]=min_cost;
        
    }
    
    isend.l=min(TSP_cost(new_chrom1,pos_s),TSP_cost(new_chrom2,pos_s));
    cout<<"Square - Rank "<<isend.rank<<" - Cost : "<<isend.l<<endl;
    MPI_Reduce(&isend, &irecv, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI::COMM_WORLD);
    if(rank==0)
        cout<<"Rank "<<irecv.rank<<" has the lowest value: Cost = "<<irecv.l<<endl;
    
    
    MPI::Finalize();

    square_pos.close();
    
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
