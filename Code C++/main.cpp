#include <iostream>
#include "matrix.cpp"
#include <vector>
#include <random>
#include <cmath>
#include <fstream>
using namespace std;

// Parametre et initialisation
const float lambda = 4.;
const float gamm = 4.;
const float lambda3[3] = {12., 12., 4.};

const int kmax = 20;
const int p = 200;
const int regimes = 3;
const float nu = 4.;

float F(float x){
    return pow(x, 4);
//    return pow(x, 10);
}

float CSA_F(float x){
    return pow(x, gamm);
}

float CSA_dF(float x){
    return pow(x, gamm-1);
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
    int k = floor(x*p);
    if (k < p) {
        return(k);
    }
    else{
        return(p);
    }

}
// x = sub/kmax*beta + (1-beta)*(sub+1)/kmax
float beta(float x, int sub){
    return (float((sub+1)-x*p));
}

inline long double Factorial(long double x) {
    if(x==0){return 1;}
    else{return (x == 1 ? x : x * Factorial(x - 1));}
}

void tirage(vector<int>& Time){
    int sum = 0;
    while (sum < p ) {
        sum += subdivision(Exponential());
        Time.push_back(sum);
    }
    Time.pop_back();
}

void sizeOrder(vector<int>& Size){
    int sum = 0;
    while (sum < p ) {
        sum += subdivision(Exponential());
        Size.push_back(sum);
    }
    Size.pop_back();
}




// a changer p et kmax
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
            float alpha = INT16_MAX;
            for(int a=1; a<k+1; a++){
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
    float alpha = F(k);
    for(int a=0; a<k+1; a++){
        if (alpha > v[k-a][T-1] + F(a)){
            alpha = v[k-a][T-1] + F(a);
        }
    }
    float h = lambda/float(p);
//    cout << h<<"  "<< alpha<< endl;
    v[k][T] = v[k][T-1] - h*(v[k][T-1]-alpha);
}

void totvect(float v[kmax][p]){
    for (int i=1; i<kmax; i++){
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
            float temp = kmax;
            float temp2 = v[kmax][0];
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


// ****** Compound Poisson Process ******
void probamu(vector<float>& mu){
    mu.push_back(exp(-nu)/(float(1-exp(-nu))));
    for (int i =1 ; i< kmax; i++){
        mu.push_back(mu[mu.size()-1]*nu/float(i));
    }
}


void PoissProcess2(float v[kmax][p], int a[kmax][p], vector<float> mu){
    
    for(int j=1; j<p; j++){
        
        // Calcul de v[i][j]
        for (int k=1; k<kmax; k++){
            float temp3 = 0;
            for (int i=1; i<a[k][j-1]+1; i++){
                temp3 += mu[i];
            }
            temp3 = 1 - temp3;
            
            float temp = (v[k-a[k][j-1]][j-1] + F(a[k][j-1]))*temp3;
            float temp2 = 0;
            for (int y=1; y<a[k][j-1]+1; y++){
                temp2 += mu[y]*(v[k-y][j-1] + F(y));
            }
            v[k][j] = v[k][j-1] + (-lambda/float(p))*(v[k][j-1] - temp - temp2);
        }
        for (int i=1; i<kmax; i++){
            float temp = kmax;
            float temp2 = F(kmax);
            for(int a=1; a<=i; a++){
                if (temp2 > v[i-a][j] + F(a)){
                    temp2 = v[i-a][j] + F(a);
                    temp = a;
                }
            }
            a[i][j] = temp;
            //            cout << temp << endl;
        }
    }
    // Calcul de a[i][p]
    for (int i=1; i<kmax; i++){
        float temp = kmax;
        float temp2 = F(kmax);
        for(int a=1; a<=i; a++){
            if (temp2 > v[i-a][p-1] + F(a)){
                temp2 = v[i-a][p-1] + F(a);
                temp = a;
            }
        }
        a[i][p-1] = temp;
    }
}








void mini(float v[kmax][p][regimes], float G[kmax][p][regimes], int i, int j){
    //assert(k>=1 && T>0);
    for (int r=0; r<regimes; r++){
        float temp = F(kmax);
        for(int a=1; a<=i; a++){
            if (temp > v[i-a][j][r] + F(a)){
                temp = v[i-a][j][r] + F(a);
            }
        }
        G[i][j][r] = v[i][j][r] - temp;
//        cout <<  G[i][j][r] << endl;
    }
}


void RSS(float v[kmax][p][regimes], float G[kmax][p][regimes], int Q[regimes][regimes]){
    for(int j=1; j<p; j++){
        for (int i=1; i<kmax; i++){
            mini(v, G, i, j-1);
        }
        // Calcul de v[i][j]
        for (int k=1; k<kmax; k++){
            for(int r=0; r<regimes; r++){
                float temp2 = 0;
                for (int y=0; y<regimes; y ++) {
                    if( y != r){
                        temp2 += Q[r][y]*(v[k][j-1][y]-v[k][j-1][r]);
//                        cout << temp2 << endl;
                    }
                }
                v[k][j][r] = v[k][j-1][r] - lambda3[r]/float(p)*G[k][j-1][r] + temp2/p;
            }
        }
    }
}




void deduireOOSRSS(float v[kmax][p][regimes], float a[kmax][p][regimes]){
    for (int r=0; r<regimes; r++){
        for (int j=0; j<p; j++){
            for (int i=1; i<kmax; i++){
                float temp = kmax;
                float temp2 = F(kmax);
                for(int a=1; a<=i; a++){
                    if (temp2 > v[i-a][j][r] + F(a)){
                        temp2 = v[i-a][j][r] + F(a);
                        temp = a;
                    }
                }
                a[i][j][r] = temp;
                cout << r << " : " << temp << endl;
            }
        }
    }
}



void write_csv(float v[kmax][p][regimes], string filename, int w){
    ofstream out(filename);
    for(int i =0; i < kmax; i++){
        for(int j = 0; j < p; j++){
            if(j <p-1){
                out << v[i][j][w]<<",";
            }
            else{
                out << v[i][j][w]<<endl;
            }
        }
    }
    out.close();
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

void write_csv(vector<int> Time, string file){
    ofstream out(file);
    for (int i = 0; i < Time.size(); i ++){
        out << Time[i] << endl;
    }
    out.close();
}

void write_csv(vector<float> Time, string file){
    ofstream out(file);
    for (int i = 0; i < Time.size(); i ++){
        out << Time[i] << endl;
    }
    out.close();
}
// ****** Simulations *******
// Montrer deux courbes : le nombre de shares restants et les ordres arrivant.

void ContinuousSaleAmount(vector<float>& a){
    a.push_back(0.5);
    for (int i=0; i<p; i++){
        float aa = a[a.size()-1];
        float temp1 = aa*(1-aa);
        temp1 = temp1*(CSA_dF(1-aa)-1);
        temp1 = lambda/float(gamm-1)*temp1;
        temp1 = aa + temp1/float(p);
        a.push_back(temp1);
    }
}



int main()
{
    // trajectoire poisson process
    InitRandom();
    vector<int> Time;
    tirage(Time);
    write_csv(Time, "TimeTEST.txt");
    

    
    // initialisation de la grille
    float v[kmax][p];
    for (int i = 0; i<kmax; i++){
        v[i][0] = F(i);
    }
    for (int i = 0; i<p; i++){
        v[0][i] = 0; // Peut etre F(0);
        v[1][i] = F(1);
    }
    
    float a[kmax][p];
    for (int i = 0; i<kmax; i++){
        a[i][0] = floor(i/2);
    }
    for (int i = 0; i<p; i++){
        a[0][i] = 0; // Peut etre F(0);
        a[1][i] = 1;
    }
    
    
    InitRandom();
    totvect(v);
    //AfficheTab(v);
    

    cout << endl;
//    cout << endl;
    deduireOOS(v, a);
    //AfficheTab(a);
//
    write_csv(v, "vTEST.txt");
    write_csv(a, "aTEST.txt");

    
    

    
    int ap[kmax][p];
    float vp[kmax][p];
    
    for (int i = 0; i<kmax; i++){
        vp[i][0] = F(i);
    }
    for (int i = 0; i<p; i++){
        vp[0][i] = 0;
    }
    
    for (int i = 0; i<kmax; i++){
        ap[i][0] = floor(i/2);
    }
    for (int i = 0; i<p; i++){
        ap[0][i] = 0; // Peut etre F(0);
        ap[1][i] = 1;
    }
    vector<float> mu;
    probamu(mu);


    PoissProcess2(vp, ap, mu);
    AfficheTab(ap);
    cout << endl;
    AfficheTab(vp);
    write_csv(vp, "vpoisTEST.txt");
    write_csv(ap, "apoisTEST.txt");
    

    
    float Vs[kmax][p][regimes];
    float G[kmax][p][regimes];
    int Q[regimes][regimes] = {{-2,1,0}, {2,-4,2}, {0,3,-2}};
    
    for (int r=0; r<regimes; r++) {
        for(int i=0; i<kmax; i++){
            Vs[i][0][r] = F(i);
        }
        for(int i=0; i<p; i++){
            Vs[0][i][r] = 0;
        }
    }
    float as[kmax][p][regimes];
    
    for (int r=0; r<regimes; r++) {
        for(int i=0; i<kmax; i++){
            as[i][0][r] = floor(i/2);
        }
        for(int i=0; i<p; i++){
            as[0][i][r] = 0;
            as[1][i][r] = 1;
        }
    }
    RSS(Vs, G, Q);
    deduireOOSRSS(Vs, as);

    
    write_csv(Vs,"Vs1TEST.txt",0);
    write_csv(Vs,"Vs2TEST.txt",1);
    write_csv(Vs,"Vs3TEST.txt",2);
    
    write_csv(as,"as1TEST.txt",0);
    write_csv(as,"as2TEST.txt",1);
    write_csv(as,"as3TEST.txt",2);

    
    
    
    
    
//********** Continuous Sale Amount *********
    vector<float> a2;
    ContinuousSaleAmount(a2);
    write_csv(a2, "acontinuousTEST.txt");
    return 0;
}



