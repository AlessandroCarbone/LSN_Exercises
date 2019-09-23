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
    
    /////////////////////////////////////////////////////////
    //CIRCUMFERENCE with R=1 and origin in (0,0)
    /////////////////////////////////////////////////////////
    
    //Test part of the code
    int N=30;       //Number of cities
    
    //Here we create a vector which will be used to create the population and to test the function of the mutations
    vector <int> order;
    for(int i=0;i<N;i++){
        order.push_back(i+1);
    }
    rnd.SaveSeed();
    
    TSP_position pos[N];
    
    //Si dispongono i punti su una circonferenza unitaria
    double angle=0;
    for(int i=0;i<N;i++){
        angle=2*M_PI*rnd.Rannyu();
        pos[i].x=cos(angle);
        pos[i].y=sin(angle);
    }
    
    ofstream circ_pos;
    circ_pos.open("Results/circ_position.txt");
    
    for(int i=0;i<N;i++){
        circ_pos<<pos[i].x<<" "<<pos[i].y<<endl;
    }
    
    //cout<<"Check: "<<TSP_check(order)<<endl;
    //cout<<"Cost: "<<TSP_cost(order,pos)<<endl;
    
    //Here we randomly create the first generation of the population
    vector <vector <int>> pop;
    int N_pop=20;       //Number of chromosomes in a population
    for(int i=0;i<N_pop;i++){
        for(int j=0;j<N;j++)
            random_shuffle(order.begin(),order.end());
        pop.push_back(order);
    }
    
    cout<<"First population"<<endl;
    for(int i=0;i<pop.size();i++){
        cout<<"Chromosome "<<i+1<<" : ";
        for(auto el: pop[i]){
            cout<<el<<" ";
        }
        cout<<endl;
    }
    cout<<endl<<endl;
    
    vector <vector <int>> new_pop;
    
    int steps=pow(10,4);
    double rand=0, index1=0, index2=0;
    
    vector <double> fitness;
    for(int i=0;i<N_pop;i++){
        fitness.push_back(0.);
    }
    
    vector <double> average_cost;  //Based on the 10 best paths
    for(int i=0;i<steps;i++){
        average_cost.push_back(0.);
    }
    
    for(int i=0;i<steps;i++){
        
        //Population fitness of the i-th generation
        for(int j=0;j<N_pop;j++){
            fitness[j]=TSP_cost(pop[j],pos);
        }
        
        //Creation of the new population
        while(new_pop.size()<N_pop){
            
            //Here we select the first and the second chromosome
            index1=selection(fitness,rnd.Rannyu());
            index2=selection(fitness,rnd.Rannyu());
            
            new_pop.push_back(pop[index1]);
            new_pop.push_back(pop[index2]);
            
            //Crossover
            if(rnd.Rannyu()<0.7)
                crossover(new_pop[new_pop.size()-1], new_pop[new_pop.size()-2], (int)round(rnd.Rannyu(1,N)));
            
            //Mutation
            //Mutation of the second chromosome
            rand=rnd.Rannyu();
            if(rand<0.05) pair_permutation(new_pop[new_pop.size()-1],(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.05 && rand<0.1)
                shift(new_pop[new_pop.size()-1],(int)round(rnd.Rannyu(1,N)));
            if(rand>0.1 && rand<0.15) invert_order(new_pop[new_pop.size()-1],(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.15 && rand<0.2)
                partial_permutation(new_pop[new_pop.size()-1], (int)round(rnd.Rannyu(0,7)), (int)round(rnd.Rannyu(7,15)), (int)round(rnd.Rannyu(15,22)), (int)round(rnd.Rannyu(22,30)));
            
            //Mutation of the first chromosome
            rand=rnd.Rannyu();
            if(rand<0.05) pair_permutation(new_pop[new_pop.size()-2],(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.05 && rand<0.1)
                shift(new_pop[new_pop.size()-2],(int)round(rnd.Rannyu(1,N)));
            if(rand>0.1 && rand<0.15) invert_order(new_pop[new_pop.size()-2],(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.15 && rand<0.2)
                partial_permutation(new_pop[new_pop.size()-2], (int)round(rnd.Rannyu(0,7)), (int)round(rnd.Rannyu(7,15)), (int)round(rnd.Rannyu(15,22)), (int)round(rnd.Rannyu(22,30)));
        }
    
        
        //Population average cost (based on the 10 best paths) of the i-th generation
        sort(fitness.begin(),fitness.end());
        for(int k=0;k<(N_pop/2);k++){
            average_cost[i]+=fitness[k];
        }
        average_cost[i]/=(N_pop/2);
        
        //Association of the new population to the current population
        //cout<<"Population "<<i+1<<endl;
        for(int k=0;k<N_pop;k++){
            pop[k]=new_pop[k];
        }
        
        //new_pop.erase(new_pop.begin(),new_pop.end());
        new_pop.clear();
    }
    
    //Here we get the best path (the chromosome of the last population which minimizes the cost)
    int best_path_index=0;
    double min_fit=TSP_cost(pop[0],pos);
    for(int i=1;i<N_pop;i++){
        if(TSP_cost(pop[i],pos)<min_fit)
            best_path_index=i;
    }
    
    circ_pos<<endl;
    for(auto el: pop[best_path_index]){
        circ_pos<<el<<endl;
    }
    circ_pos.close();
    
    //Data Blocking
    vector <double> sum_prog;
    vector <double> err_prog;
    
    for(int l=0;l<steps;l++){
        sum_prog.push_back(0);
        err_prog.push_back(0);
    }
    
    data_blocking(average_cost,sum_prog,err_prog);
    
    ofstream circ_output;
    circ_output.open("Results/circ_output.txt");
    
    for(int i=0;i<steps;i++){
        circ_output<<average_cost[i]<<" "<<sum_prog[i]<<" "<<err_prog[i]<<endl;
    }
    
    circ_output.close();
    
    rnd.SaveSeed();
    
    /////////////////////////////////////////////////////////
    //SQUARE with L=1, x in (0,1), y in (0,1)
    /////////////////////////////////////////////////////////
    
    
    for(int i=0;i<N;i++){
        rand=rnd.Rannyu();
        pos[i].x=rand;
        rand=rnd.Rannyu();
        pos[i].y=rand;
    }
    
    ofstream square_pos;
    square_pos.open("Results/square_position.txt");
    
    for(int i=0;i<N;i++){
        square_pos<<pos[i].x<<" "<<pos[i].y<<endl;
    }
    
    for(int i=0;i<N_pop;i++){
        for(int j=0;j<N;j++)
            random_shuffle(order.begin(),order.end());
        pop[i]=order;
    }
    
    
    for(int i=0;i<steps;i++){
        
        //Population fitness of the i-th generation
        for(int j=0;j<N_pop;j++){
            fitness[j]=TSP_cost(pop[j],pos);
        }
        
        //Creation of the new population
        while(new_pop.size()<N_pop){
            
            //Here we select the first and the second chromosome
            index1=selection(fitness,rnd.Rannyu());
            index2=selection(fitness,rnd.Rannyu());
            
            new_pop.push_back(pop[index1]);
            new_pop.push_back(pop[index2]);
            
            //Crossover
            if(rnd.Rannyu()<0.7)
                crossover(new_pop[new_pop.size()-1], new_pop[new_pop.size()-2], (int)round(rnd.Rannyu(1,N)));
            
            //Mutation
            //Mutation of the second chromosome
            rand=rnd.Rannyu();
            if(rand<0.05) pair_permutation(new_pop[new_pop.size()-1],(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.05 && rand<0.1)
                shift(new_pop[new_pop.size()-1],(int)round(rnd.Rannyu(1,N)));
            if(rand>0.1 && rand<0.15) invert_order(new_pop[new_pop.size()-1],(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.15 && rand<0.2)
                partial_permutation(new_pop[new_pop.size()-1], (int)round(rnd.Rannyu(0,7)), (int)round(rnd.Rannyu(7,15)), (int)round(rnd.Rannyu(15,22)), (int)round(rnd.Rannyu(22,30)));
            
            //Mutation of the first chromosome
            rand=rnd.Rannyu();
            if(rand<0.05) pair_permutation(new_pop[new_pop.size()-2],(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.05 && rand<0.1)
                shift(new_pop[new_pop.size()-2],(int)round(rnd.Rannyu(1,N)));
            if(rand>0.1 && rand<0.15) invert_order(new_pop[new_pop.size()-2],(int)round(rnd.Rannyu(0,N-1)),(int)round(rnd.Rannyu(0,N-1)));
            if(rand>0.15 && rand<0.2)
                partial_permutation(new_pop[new_pop.size()-2], (int)round(rnd.Rannyu(0,7)), (int)round(rnd.Rannyu(7,15)), (int)round(rnd.Rannyu(15,22)), (int)round(rnd.Rannyu(22,30)));
        }
        
        
        //Population average cost (based on the 10 best paths) of the i-th generation
        sort(fitness.begin(),fitness.end());
        for(int k=0;k<(N_pop/2);k++){
            average_cost[i]+=fitness[k];
        }
        average_cost[i]/=(N_pop/2);
        //cout<<average_cost[i]<<endl;
        
        //Association of the new population to the current population
        //cout<<"Population "<<i+1<<endl;
        for(int k=0;k<N_pop;k++){
            pop[k]=new_pop[k];
        }
        
        //new_pop.erase(new_pop.begin(),new_pop.end());
        new_pop.clear();
    }
    
    //Here we get the best path (the chromosome of the last population which minimizes the cost)
    best_path_index=0;
    min_fit=TSP_cost(pop[0],pos);
    for(int i=1;i<N_pop;i++){
        if(TSP_cost(pop[i],pos)<min_fit)
            best_path_index=i;
    }
    
    square_pos<<endl;
    for(auto el: pop[best_path_index]){
        square_pos<<el<<endl;
    }
    square_pos.close();
    
    //Data Blocking
    data_blocking(average_cost,sum_prog,err_prog);
    
    ofstream square_output;
    square_output.open("Results/square_output.txt");
    
    for(int i=0;i<steps;i++){
        square_output<<average_cost[i]<<" "<<sum_prog[i]<<" "<<err_prog[i]<<endl;
    }
    
    square_output.close();
    
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
