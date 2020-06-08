hmm3.cpp
Type
C++
Size
6 KB (5,639 bytes)
Storage used
0 bytes
Location
hmm3 final
Owner
Theodor Panagiotakopoulos
Modified
Sep 12, 2019 by Theodor Panagiotakopoulos
Created
Sep 12, 2019
Add a description
Viewers can download

#include <cstdlib>
#include <iostream>
#include <string.h>
#include <tuple>
#include<cmath>
#include "vec.h"
using namespace std;


double* getinput(bool f = true){
    double rows, cols;
    if(f){
        cin >> rows;
    }else{
        rows = 1;
    }
    cin >> cols;
    double *out = new double[(int)(rows*cols+2)];
    out[0] = rows;
    out[1] = cols;
    for (int i=2;i<rows*cols+2;i++){
        cin >> out[i];
    }
    return out;
}


vec a_pass(vec &A,vec &B, vec &P, vec &obs){
    int ns = A.rows;
    vec ALPHA (obs.cols,ns);
    ALPHA.insert_row(P.trans()^(B.get_col(obs(0,0))),0);
    //vec tmp = elmnt_wise(P,extrct_col(B,ob.pi[0]));
    //copy_line(ALPHA,tmp,0,0);
    
    for(int i=1; i<obs.cols; i++){
        ALPHA.insert_row((ALPHA.get_row(i-1)*A).trans()^B.get_col(obs(0,i)),i);
        //vec tmp1 = mat_mul(extrct_row(ALPHA,i-1),A);
        //copy_line(ALPHA,elmnt_wise(tmp1,extrct_col(B,ob.pi[i])),i,0);
    }
    return ALPHA;
}
pair <vec, vec> a_pass_scaled(vec &A,vec &B, vec &P, vec &obs){
    int ns = A.rows;
    vec ALPHA (obs.cols,ns);
    
    vec tmp = P.trans()^(B.get_col(obs(0,0)));
    vec c (obs.cols,1);
    // find c0;
    for(int i=0;i<ns;i++){
        c(0,0)+=tmp(0,i);
    }
    c(0,0)=1/c(0,0);
    //scale ALPHA0
    for(int i=0;i<ns;i++){
        ALPHA(0,i)=c(0,0)*tmp(0,i);
    }
    for(int i=1; i<obs.cols; i++){
        vec tmp1 = (ALPHA.get_row(i-1)*A)^(B.get_col(obs(0,i)).trans());
        for(int j=0;j<ns;j++){
            c(i,0)+=tmp1(0,j);
        }
        c(i,0)=1.0/c(i,0);
        for(int j=0;j<ns;j++){
            ALPHA(i,j)=c(i,0)*tmp1(0,j);
        }
    }
    pair <vec, vec> out (ALPHA,c);
    return out;
}
vec b_pass(vec & A, vec & B, vec & obsv){
    int T = obsv.cols;
    vec OUT_betas (T - 1, A.rows); // T must be greater than 2
    for(int i=0;i<A.rows;i++){
        OUT_betas(T-2,i)=1;
    }
    for(int i=T-3;i>=0;i--){
        OUT_betas.insert_row((OUT_betas.get_row(i+1)^(B.get_col(obsv(0,i+2)).trans()))*A.trans(),i);
    }
    return OUT_betas;
}

vec b_pass_scaled(vec & A, vec & B, vec & obsv, vec & c){
    int T = obsv.cols;
    vec OUT_betas (T - 1, A.rows); // T must be greater than 2
    for(int i=0;i<A.rows;i++){
        OUT_betas(T-2,i)=c(T-1,0);
    }
    for(int i=T-3;i>=0;i--){
        vec tmp = (OUT_betas.get_row(i+1)^(B.get_col(obsv(0,i+2)).trans()))*A.trans();
        for(int j=0;j<A.rows;j++){
            OUT_betas(i,j)=tmp(0,j)*c(i+1,0);
        }
    }
    return OUT_betas;
}
pair <vec, vec> gamma(vec & A, vec & B, vec & p,vec & obsv, vec &alpha, vec & beta){
    double sum=0;
    for(int i=0;i<A.rows;i++){
        sum+=alpha(alpha.rows-1,i);
    }
    vec gamma_sum (A.rows,A.rows);
    vec gammai (obsv.cols-1, A.rows);
    vec tmp (A.rows,A.cols);
    for (int i=0;i<obsv.cols-1;i++){
        tmp = ((alpha.get_row(i).trans()*(beta.get_row(i)^(B.get_col(obsv(0,i+1)).trans())))^A);
        tmp/=sum;
        gammai.insert_row(tmp.sum_cols().trans(),i);
        gamma_sum+=tmp;
    }
    pair<vec, vec> out (gammai, gamma_sum);
    return out;
    
}
tuple <vec, vec, vec, vec> estimation(vec & A, vec & B, vec & p,vec & obsv){ //estimation of A , B, pi
    //calculate alpha beta
    pair <vec, vec> outalpha= a_pass_scaled(A,B,p,obsv);
    vec alpha = outalpha.first;
    vec c = outalpha.second;
    vec beta = b_pass_scaled(A,B,obsv,c);

    pair <vec, vec> out = gamma(A, B, p, obsv,alpha, beta);
    vec gammai = out.first;
    vec digamma_sum = out.second;
    vec newA (A.rows,A.cols);
    for(int i=0;i<A.rows;i++){
        double rowsum = 0;
        for(int j=0;j<A.rows;j++){
            rowsum+= digamma_sum(i,j);
        }
        for(int j=0;j<A.rows;j++){
            newA(i,j)= digamma_sum(i,j)/rowsum;
        }
    }
    vec newB(B.rows,B.cols);
    vec obs_over_t (obsv.cols -1,B.cols);
    for(int i=0;i<obsv.cols-1;i++){
        obs_over_t(i,obsv(0,i)) = 1;
    }
    vec Bdraft = gammai.trans()*obs_over_t;
    vec tmpsum(1,obsv.cols-1);
    for(int i=0;i<obsv.cols-1;i++){
        tmpsum+=gammai.get_row(i);
    }
    for(int i=0;i<B.rows;i++){
        for(int j=0;j<B.cols;j++){
            newB(i,j) = Bdraft(i,j)/tmpsum(0,i);
        }
    }
    vec newpi = gammai.get_row(0);

    tuple<vec, vec, vec, vec> out2 (newA,newB,newpi, c);

    return out2;
}
vec fix_vec(vec & A){
    vec B (A.rows,A.cols);
    for(int i=0;i<A.rows;i++){
        for(int j=0;j<A.cols;j++){
            if (A(i,j) < 0.0001){
                B(i,j)=0;
            }else{
                B(i,j)=A(i,j);
            }
        }
    }
    return B;
}

tuple <vec,vec,vec> baum_welch(vec & A, vec & B, vec & p,vec & obsv,int maxiter){
    tuple <vec, vec, vec, vec> out = estimation(A, B, p, obsv);
    vec newA = get<0>(out);
    vec newB = get<1>(out);
    vec newp = get<2>(out);
    int i=0;
    double logprob = -1;
    double oldprob = -1;
    while(i<maxiter&& (logprob > oldprob || i<2)){
        oldprob = logprob;
        out = estimation(newA, newB, newp, obsv);
        vec tmp = get<0>(out);
        newA = fix_vec(get<0>(out));
        newB = fix_vec(get<1>(out));
        newp = fix_vec(get<2>(out));
        vec c = get<3>(out);

        logprob=0;
        for(int j=0;j<c.rows;j++){
            logprob+=log10(c(j,0));
        }
        logprob=-logprob;
        
        i++;
    }
    tuple <vec,vec,vec> o (newA,newB,newp);
    return o;
    
}
int main(){
    vec A (getinput());
    vec B (getinput());
    vec p (getinput());
    vec obsv (getinput(false));
    tuple <vec,vec,vec> baum = baum_welch(A, B, p,obsv,1040);
    get<0>(baum).print_inline();
    get<1>(baum).print_inline();
    return 0;
}
