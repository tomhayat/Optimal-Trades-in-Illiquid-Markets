#include <iostream>
#include "matrix.h"
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
using namespace std;
const float lambda = 1.;
const int kmax= 200;
const int p = 1000;
const float nu = 1;
float F(float x){
    return x*x/2;
}

// Simulation de lois
void InitRandom(){
    srand((unsigned int)time(0));
}
float Uniform(){
    return (float(rand())/(float(RAND_MAX)+1.));
}
float Exponential(){
    return (-log(Uniform())/lambda);
}

// reperer un reel dans la subdivision
int subdivision(float x){
    int k = floor(x*kmax);
    return(k);
}
// x = sub/kmax*beta + (1-beta)*(sub+1)/kmax
float beta(float x, int sub){
    return (float((sub+1)-x*kmax));
}

// renvoie v_k_t en ayant une partie de la grille
void v_k_t(float v[p][kmax], int k, int T){
    assert(k > 1 && T > 0);
    float vkT = 0.;
    // Methode de Monte Carlo
    int iteration = 100000; // optimal ?
    for (int i=0; i <iteration; i++){
        // On tire le sigma au hasard
        float sigma_simul = Exponential();
        int sigma = subdivision(sigma_simul);
        float tau = beta(sigma_simul,sigma);
        // sigma_simul*kmax = sigma*tau + (1-tau)*(sigma+1)
        float temp = 0.;
        // On regarde ces 2 cas pour sigma et (sigma +1)
        
        if (sigma >= T) {
            temp += F(k)*tau;
        }
        else{
            float alpha= INT16_MAX;
            for(int a=1; a < k+1; a++){
                if (alpha > v[T-sigma][k-a] + F(a)){
                    alpha = v[T-sigma][k-a] + F(a);
                }
            }
            temp += alpha*tau;
        }
        if ((sigma+1) >= T) {
            temp += F(k)*(1-tau);
        }
        else{
            float alpha= INT16_MAX;
            for(int a=1; a < k+1; a++){
                if (alpha > v[T-(sigma+1)][k-a] + F(a)){
                    alpha = v[T-(sigma+1)][k-a] + F(a);
                }
            }
            temp += alpha*(1-tau);
        }
        vkT += temp;
    }
    v[T][k] = vkT/iteration;
}



void v_k_t_diff(float v[kmax][p], int k, int T){
    //assert(k>=1 && T>0);
    float alpha = INT16_MAX;
    for(int a=0; a<k+1; a++){
        if (alpha > v[k-a][T-1] + F(a)){
            alpha = v[k-a][T-1] + F(a);
        }
    }
    v[k][T] = v[k][T-1] - lambda/p*(v[k][T-1]-alpha);
}

void totvect(float v[kmax][p]){
    for (int i=2; i<kmax; i++){
        for (int j=1; j<p; j++){
            v_k_t_diff(v, i, j);
        }
    }
}



void AfficheTab(float v[kmax][p]){
    for (int i=0; i<kmax; i++){
        for (int j=0; j<p; j++){
            cout << v[i][j] << "  ";
        }
        cout << endl;
    }
}

void AfficheTab(int v[kmax][p]){
    for (int i=0; i<kmax; i++){
        for (int j=0; j<p; j++){
            cout << v[i][j] << "  ";
        }
        cout << endl;
    }
}

void deduireOOS(float v[kmax][p], float a[kmax][p]){
    
    for (int j=0; j<p; j++){
        for (int i=1; i<kmax; i++){
            float temp = INT16_MAX;
            float temp2 = INT16_MAX;
            for(int a=1; a<=i; a++){
                if (temp2 > v[i-a][j] + F(a)){
                    temp2 = v[i-a][j] + F(a);
                    temp = a;
                }
            }
            a[i][j] = temp;
        }
    }
}

void write_csv(float v[kmax][p], string file){
    ofstream out(file);
    for(int i =0; i < kmax; i++){
        for(int j = 0; j < p; j++){
            if(j <p-1){
            out << v[i][j]<<",";
            }
            else{
                out << v[i][j]<<endl;
            }
        }
    }
    out.close();
}

void write_csv(int v[kmax][p], string file){
    ofstream out(file);
    for(int i =0; i < kmax; i++){
        for(int j = 0; j < p; j++){
            if(j <p-1){
                out << v[i][j]<<",";
            }
            else{
                out << v[i][j]<<endl;
            }
        }
    }
    out.close();
}

// ****** Compound Poisson Process ******
// Ne pas oublier que cest v~.
void PoissProcess(float v[kmax][p], int a[kmax][p]){
    
    for(int j =1; j <p; j++){
        for (int i=1; i<kmax; i++){
            float temp = INT16_MAX;
            float temp2 = INT16_MAX;
            for(int a=1; a<=i; a++){
                if (temp2 > v[i-a][j-1] + F(a)){
                    temp2 = v[i-a][j-1] + F(a);
                    temp = a;
                }
            }
            a[i][j-1] = temp;
        }
        // Calcul de v[i][j]
        for (int k = 1; k < kmax; k ++){
            float temp = (v[k-a[k][j-1]][j-1] + F(a[k][j-1]))*exp(-nu*a[k][j-1]);
            float temp2 = 0;
            for (int y = 1; y < a[k][j-1]+1; y ++) {
                temp2 += nu*exp(-nu*y)*(v[k-y][j-1] + F(y));
            }
            v[k][j] = v[k][j-1] - lambda/p*(v[k][j-1] - temp + temp2);
        }
    }
    // Calcul de a[i][p]
    for (int i=1; i<kmax; i ++){
        float temp = INT16_MAX;
        float temp2 = INT16_MAX;
        for(int a = 1; a <= i; a ++){
            if (temp2 > v[i-a][p-1] + F(a)){
                temp2 = v[i-a][p-1] + F(a);
                temp = a;
            }
        }
        a[i][p-1] = temp;
    }
    
    
}


int main()
{
    // initialisation de la grille
    float v[kmax][p];
    for (int i = 0; i<kmax; i++){
        v[i][0] = F(i);
    }
    for (int i = 0; i<p; i++){
        v[0][i] = 0; // Peut etre F(0);
        v[1][i] = F(1);
    }
    
    
    InitRandom();
    totvect(v);
    //AfficheTab(v);
    
    
    float a[kmax][p];
    //    for (int i = 0; i<kmax; i++){
    //        a[i][0] = F(i);
    //    }
    //    for (int i = 0; i<p; i++){
    //        a[0][i] = 0; // Peut etre F(0);
    //        a[1][i] = F(1);
    //    }
    cout << endl;
    cout << endl;
    deduireOOS(v, a);
    //AfficheTab(a);
    write_csv(v, "v.txt");
    write_csv(a,"a.txt");
    
    float vc[kmax][p];
    for (int i = 0; i<kmax; i++){
        vc[i][0] = F(i);
    }
    for (int i = 0; i<p; i++){
        vc[0][i] = 0;
    }
    int ac[kmax][p];
    PoissProcess(vc, ac);
    AfficheTab(vc);
    AfficheTab(ac);
    
    write_csv(vc, "vpoisscomp.txt");
    write_csv(ac, "apoisscomp.txt");
    return 0;
}



