hmm3_v1.cpp
Type
C++
Size
18 KB (18,305 bytes)
Storage used
0 bytes
Location
AI_A1
Owner
Theodor Panagiotakopoulos
Modified
Sep 11, 2019 by Theodor Panagiotakopoulos
Opened
11:42 PM by me
Created
Sep 10, 2019
Add a description
Viewers can download

/////////////////////
// Assignment 1 - AI
////////////////////

#include <cstdlib>
#include <iostream>
#include <string.h>
#include <tuple>
#include<cmath>
using namespace std;

class vec{ //create a new type of variable called vec
    public:
        int rows;
        int cols;
        double **pi;
        vec(){

        }
        ~vec(){
            if(rows>0 && cols>0){
                for(int i=0;i<rows;i++){
                    //cout<< rows;
                    delete[] pi[i];
                }
                delete[] pi;
            }
        }
        vec(const vec& other){
            operator=(other);            
        }
        vec& operator=(const vec &other){
            rows = other.rows;
            cols = other.cols;
            pi = new double*[other.rows];
            for(int i=0;i<rows;i++){
                pi[i] = new double[other.cols];
                for(int j=0;j<other.cols;j++){
                   pi[i][j] = other.pi[i][j];
                }
                //memcpy( &pi[i], &other.pi[i], sizeof(other.pi[i]));
            }
            return *this;
        }
        void initialize(int x, int y, double val=0){
            rows = x;
            cols = y;
            pi = new double*[x];
            for(int i=0;i<x;i++){
                pi[i] = new double[y];
                for(int j=0;j<y;j++){
                    pi[i][j]=val;
                }
            }
        }
        double*& operator[](int idx){
            return pi[idx];
        }
        vec operator/(double x){
            for(int i=0;i<rows;i++){
                for(int j=0;j<cols;j++){
                    pi[i][j] /=x;
                }
            }
            return *this;
        }
        vec operator +=(vec B){
            for(int i=0;i<rows;i++){
                for(int j=0;j<cols;j++){
                    pi[i][j] +=B[i][j];
                }
            }
            return *this;
        }
        /*double& operator[](int idx, int idx2){
            return pi[idx][idx2];
        }*/
};

class multi_vec{
    public:
        int x1;
        int x2;
        int x3;
        vec *l;
        multi_vec(int t,int n, int m){
            l = new vec[t];
            x1=t,x2=n;x3=m;
            for(int i=0;i<t;i++){
                l[i].initialize(n,m);
            }
        }
        vec& operator[](int idx){
            return l[idx];
        }
};


////////
// HMM0
////////

vec mat_mul(vec v_in1, vec v_in2){//input 1 is transition matrix and input 2 is emission matrix
    //if(v_in1.cols == v_in1.rows){
    //tratar excecao    
        vec v_out; //create vector-object output
        v_out.pi = new double *[v_in1.rows]; //adress pi of vector-output to ohave the same number of slots as the number of rows of input 1 
        v_out.rows = v_in1.rows;
        v_out.cols = v_in2.cols; //remember: A(axb)*B(bxc)=C(axc)

        for(int i=0; i<v_in1.rows; i++){ //

            v_out.pi[i] = new double[v_in2.cols]; //
            for(int j=0; j<v_in2.cols; j++){ //

                double lin_comb = 0; //
                for (int k=0;k<v_in1.cols;k++){ //

                    lin_comb+= v_in1.pi[i][k]*v_in2.pi[k][j]; //
                }
                v_out.pi[i][j] = lin_comb; //
            }
        }
        return v_out;
    //}
}

void print_array_inline(vec A){
    cout << A.rows << ' ' << A.cols;
    for (int i=0; i<A.rows; i++)    {
        for(int j=0; j< A.cols; j++){
            cout << ' ' << A.pi[i][j];
        }
    }
}

void print_vec(vec A){
    for(int i=0; i<A.rows; i++){
        for(int j=0; j<A.cols; j++){
            cout << A.pi[i][j] << " ";
        }
        cout << endl;
    }
}

vec getinput(){
    string line;
    double input[5000]={};
    int start = 0;
    int counter = 0;
    
    getline(cin, line);
    line += " ";
    bool flag = true;

    for (int i=0; i<line.length(); i++){
        if (line.at(i) == ' ' && flag){
            input[counter] = stof(line.substr(start,i-start));
            start = i+1;
            counter++;
            flag = false;
        }else if (line.at(i) != ' '){
            flag = true;
        }
    }
    vec A;
    A.rows = (int)input[0];
    A.cols = (int)input[1];
    A.pi = new double *[A.rows];
    for(int i=0; i<A.rows; i++){
        A.pi[i] = new double[A.cols];
        for(int j=0; j<A.cols; j++){
            A.pi[i][j] = input[2 + i*A.cols+ j];
        }
    }
    return A;   
}

////////
// HMM1
////////

struct list{
    int *pi;
    int size;
};

list get_observations(){
    int nums;
    cin >> nums;
    list o;
    o.size = nums;
    o.pi = new int[nums];
    for (int i=0;i<nums;i++){
        cin >> o.pi[i];
    }
    return o;
}

void print_list(list a){
    cout << '[';
    for(int i=0;i<a.size - 1;i++){
        cout << a.pi[i] << ' ';
    }
    if(a.size>0){cout << a.pi[a.size-1];}
    cout << ']';
}

vec elmnt_wise(vec in1,vec in2){
    //note: in1 and in2 must have the same dimension -> lack exception implementation//
    int n = in1.rows;
    int m = in1.cols;
    vec out; 
    out.pi = new double *[in1.rows]; 
    out.rows = n;
    out.cols = m; 

    for(int i=0; i<n; i++){
        out.pi[i] = new double[m]; //
        for(int j=0; j<m; j++){
            out.pi[i][j] = in1.pi[i][j]*in2.pi[j][i]; //remember: element wise is row by column 
        }
    }
    return out;
}

vec extrct_col(vec A, int loc){
    vec B;
    B.rows = A.rows;
    B.cols = 1;
    B.pi = new double *[A.rows];
    for(int i=0;i<A.rows;i++){
        B.pi[i] = new double[B.cols];
        B.pi[i][0] = A.pi[i][loc];
    }
    return B;

}
void copy_line(vec &A, vec B,int a, int b){
    for(int i=0;i<A.cols;i++){
        A[a][i] = B[b][i];
    }
}

vec extrct_row(vec A, int loc){
    vec B;
    B.rows = 1;
    B.cols = A.cols;
    B.pi = new double *[1];
    B.pi[0] = new double[B.cols];
    for(int i=0;i<A.cols;i++){
        B.pi[0][i] = A.pi[loc][i];
    }
    return B;
}

vec a_pass(list ob,vec A,vec B,vec P){
    int ns = A.rows;
    vec ALPHA;
    ALPHA.initialize(ob.size,ns);
    vec tmp = elmnt_wise(P,extrct_col(B,ob.pi[0]));
    copy_line(ALPHA,tmp,0,0);
    for(int i=1; i<ob.size; i++){
        vec tmp1 = mat_mul(extrct_row(ALPHA,i-1),A);
        copy_line(ALPHA,elmnt_wise(tmp1,extrct_col(B,ob.pi[i])),i,0);
    }
    return ALPHA;
}

pair <vec, vec> a_pass_norm(list ob,vec A,vec B,vec P){
    int ns = A.rows;
    vec ALPHA;
    ALPHA.initialize(ob.size,ns);
    vec tmp = elmnt_wise(P,extrct_col(B,ob.pi[0]));
    vec c;
    // find c0;
    c.initialize(ob.size,1);
    for(int i=0;i<ns;i++){
        c[0][0]+=tmp[0][i];
    }
    c[0][0]=1/c[0][0];
    //scale ALPHA0
    for(int i=0;i<ns;i++){
        ALPHA[0][i]=c[0][0]*tmp[0][i];
    }
    for(int i=1; i<ob.size; i++){
        vec tmp1 = mat_mul(extrct_row(ALPHA,i-1),A);
        tmp1 = elmnt_wise(tmp1,extrct_col(B,ob.pi[i]));
        for(int j=0;j<ns;j++){
            c[i][0]+=tmp1[0][j];
        }
        c[i][0]=1.0/c[i][0];
        for(int j=0;j<ns;j++){
            ALPHA[i][j]=c[i][0]*tmp1[0][j];
        }
    }
    pair <vec, vec> out (ALPHA,c);
    return out;
}

////////
// HMM2
////////
vec inv(vec A){
    vec B;
    B.cols = A.rows;
    B.rows = A.cols;
    B.pi = new double*[A.cols];
    for(int j=0;j<A.cols;j++){
        B.pi[j] = new double[A.rows];
        for(int i=0;i<A.rows;i++){
            B.pi[j][i]=A.pi[i][j];
        }
    }
    return B;
}
pair <vec, vec> maximus(vec input){
    //pair <vec, vec> OUT;
    vec out_pr,out_indx;
    out_pr.rows = input.rows;
    out_pr.cols = 1;
    out_pr.pi = new double *[input.rows];
    out_indx.rows = input.rows;
    out_indx.cols = 1;
    out_indx.pi = new double *[input.rows];
    int bin;
    for(int i=0;i<input.rows;i++){
        out_pr.pi[i] = new double[1];
        out_indx.pi[i] = new double[1];
        out_pr.pi[i][0] = input.pi[i][0];
        bin= 1;
        out_indx.pi[i][0] = 1;
        for(int j=1;j<input.cols;j++){
            bin*=2;
            if(out_pr.pi[i][0] < input.pi[i][j]){
               out_pr.pi[i][0] = input.pi[i][j];
               out_indx.pi[i][0] = bin;
            }else if(out_pr.pi[i][0] == input.pi[i][j]){
               out_indx.pi[i][0] += bin;
            }
        }
    }
    pair <vec, vec> OUT (out_pr,out_indx);
    return OUT;
}
void addrows(vec &A, vec B){
    double **newpi = new double *[A.rows + B.rows];
    for(int i=0;i< A.rows;i++){
        newpi[i] = A.pi[i];
    }
    for (int i=0;i<B.rows;i++){
        newpi[A.rows+i] = new double[B.cols];
        for (int j=0;j<B.cols;j++){
            newpi[A.rows + i][j] = B.pi[i][j];
        }
        
    }
    A.rows += B.rows;
    A.pi = newpi;
}
vec find_prev_state(vec di,int time,int start){
    vec out;
    out.rows = 0;
    out.cols = time;
    //out.pi[0] = new double[di.rows];
    for(int i=0; start>0; i++){
        if(start%2){
            vec tmp;
            if (time == 1){
                tmp.rows=1; tmp.cols=1;tmp.pi = new double *[1];tmp.pi[0] = new double[1];
                tmp.pi[0][0] = i + 0.0;
            }else{
                vec prev = find_prev_state(di, time-1, di.pi[time-1][i]);
                tmp.rows=prev.rows; tmp.cols=time;tmp.pi = new double *[tmp.rows]; // initialization
                for(int ik=0;ik<tmp.rows;ik++){
                    tmp.pi[ik] = new double[time];
                    tmp.pi[ik][tmp.cols-1] = i + 0.0;
                    for(int j=0;j<tmp.cols-1;j++){
                        tmp.pi[ik][j] = prev.pi[ik][j];
                    }
                }
            }
            //print_vec(tmp);
            addrows(out,tmp);
        }   
        start/=2;  
    }
    return out;    
}


vec Viterbi(vec A, vec B, vec p, list obsv){
    
    vec d = elmnt_wise(p,extrct_col(B,obsv.pi[0]));
    //print_vec(d);
    //cout << endl;
    vec dindx = extrct_row(p,0);
    for(int i=1;i<obsv.size;i++){
        /*
        cout << "--------" << endl;
        print_vec(elmnt_wise(mat_mul(extrct_col(B,obsv.pi[i]), extrct_row(d,i-1)),A));
        cout << "||" << endl;;
        print_vec(extrct_row(d,i-1));
        cout << "||" << endl;;
        print_vec(extrct_col(B,obsv.pi[i]));
        cout << "||" << endl;;
        print_vec(mat_mul(extrct_col(B,obsv.pi[i]), extrct_row(d,i-1)));
        cout << "--------" << endl;
        */
        pair<vec,vec> ans = maximus(elmnt_wise(mat_mul(extrct_col(B,obsv.pi[i]), extrct_row(d,i-1)),A));
        //cout << "--------" << endl;
        //print_vec(ans.first);
        //cout << "--------" << endl;
        addrows(d,inv(ans.first));
        addrows(dindx,inv(ans.second));
    }
    //print_vec(d);
    //print_vec(dindx);

    pair<vec,vec> ans = maximus(extrct_row(d,d.rows-1));
    
    return find_prev_state(dindx, dindx.rows, ans.second.pi[0][0]);
    //print_vec(ans.second);
    return d;

}

////////
// HMM3
////////
vec sum_cols(vec A){
    vec B;
    B.initialize(A.rows, 1,0);
    for(int i=0;i<A.rows;i++){
        for(int j=0;j<A.cols;j++){
            B[i][0] += A[i][j];
        }
        
    }
    return B;
    

}
vec b_pass(vec A, vec B, list obsv, int T){
    
    vec OUT_betas;
    OUT_betas.initialize(T - 2, A.rows); // T must be greater than 2
    vec last_line;
    last_line.initialize(1, A.rows, 1);
    addrows(OUT_betas,last_line); // this is time T
    for(int i=T-3;i>=0;i--){
        //OUT_betas.pi[i] = elmnt_wise(extrct_row(OUT_betas,i+1),mat_mul(A,extrct_col(B, obsv.pi[i+1])));
        copy_line(OUT_betas,mat_mul(elmnt_wise(extrct_row(OUT_betas,i+1),extrct_col(B, obsv.pi[i+2])), inv(A)),i,0);
    }
    return OUT_betas;
}
vec b_pass_norm(vec A, vec B, list obsv, int T, vec c){
    
    vec OUT_betas;
    OUT_betas.initialize(T - 2, A.rows); // T must be greater than 2
    vec last_line;
    last_line.initialize(1, A.rows, c[T-1][0]);
    addrows(OUT_betas,last_line); // this is time T
    for(int i=T-3;i>=0;i--){
        //OUT_betas.pi[i] = elmnt_wise(extrct_row(OUT_betas,i+1),mat_mul(A,extrct_col(B, obsv.pi[i+1])));
        vec tmp = mat_mul(elmnt_wise(extrct_row(OUT_betas,i+1),extrct_col(B, obsv.pi[i+2])), inv(A));
        for(int j=0;j<A.rows;j++){
            OUT_betas[i][j]=tmp[0][j]*c[i][0];
        }
    }
    return OUT_betas;
}

tuple<multi_vec, vec, vec> gamma(vec A, vec B, vec p,list obsv){
    //important multi_vec digamma (obsv.size-1,A.rows,A.rows);
    multi_vec digamma (0,0,0);
    pair <vec, vec> outalpha= a_pass_norm(obsv,A,B,p);
    vec alpha = outalpha.first;
    vec c = outalpha.second;
    vec beta = b_pass_norm(A,B,obsv, obsv.size,c);
    //print_vec(beta);
    //cout<< "----end----" << endl;
    //print_vec(c);
    //cout<< "----end----" << endl;
    //print_vec(alpha);
    double sum=0;
    for(int i=0;i<A.rows;i++){
        sum+=alpha[alpha.rows-1][i];
    }
    vec gamma_sum;
    vec gammai;
    gammai.initialize(0, A.rows);
    gamma_sum.initialize(A.rows,A.rows);
    vec tmp;
    tmp.initialize(A.rows,A.cols);
    for (int i=0;i<obsv.size-1;i++){
        //important digamma[i] = elmnt_wise(mat_mul(inv(extrct_row(alpha,i)),elmnt_wise(extrct_row(beta,i),extrct_col(B,obsv.pi[i+1]))),inv(A))/sum;
        tmp = elmnt_wise(mat_mul(inv(extrct_row(alpha,i)),elmnt_wise(extrct_row(beta,i),extrct_col(B,obsv.pi[i+1]))),inv(A))/sum;
        //cout << "-----digamma----" << endl;
        //print_vec(digamma[i]);
        //cout << "---------" << endl;
        addrows(gammai,inv(sum_cols(tmp)));
        gamma_sum+=tmp;
        //print_vec(gammai);
        //cout << "---------" << endl;
    }
    tuple<multi_vec, vec, vec> out (digamma,gammai, gamma_sum);
    return out;
    
}

tuple <vec, vec, vec> estimation(vec A, vec B, vec p,list obsv){ //estimation of A , B, pi
    tuple <multi_vec, vec, vec> out = gamma(A, B, p, obsv);
    multi_vec digamma = get<0>(out);
    vec gammai = get<1>(out);
    vec digamma_sum = get<2>(out);
    vec newA;
    newA.initialize(A.rows,A.cols);
    for(int i=0;i<A.rows;i++){
        double rowsum = 0;
        for(int j=0;j<A.rows;j++){
            rowsum+= digamma_sum[i][j];
        }
        for(int j=0;j<A.rows;j++){
            newA[i][j]= digamma_sum[i][j]/rowsum;
        }
    }
    vec newB;
    newB.initialize(B.rows,B.cols);
    vec obs_over_t;
    obs_over_t.initialize(obsv.size -1,B.cols);
    for(int i=0;i<obsv.size-1;i++){
        obs_over_t[i][obsv.pi[i]] = 1;
    }
    vec Bdraft;
    Bdraft = mat_mul(inv(gammai),obs_over_t);
    vec tmpsum;
    tmpsum.initialize(1,obsv.size-1);
    for(int i=0;i<obsv.size-1;i++){
        tmpsum+=extrct_row(gammai,i);
    }
    for(int i=0;i<B.rows;i++){
        for(int j=0;j<B.cols;j++){
            newB[i][j] = Bdraft[i][j]/tmpsum[0][i];
        }
    }
    vec newpi;
    newpi = extrct_row(gammai,0);
    tuple<vec, vec, vec> out2 (newA,newB,newpi);

    return out2;
}
double diff(vec A, vec B){
    double sum=0;
    for(int i=0;i<A.rows;i++){
        for(int j=0;j<A.cols;j++){
            sum += fabs((A[i][j]-B[i][j]));
        }
    }
    return sum;
}
vec fix_vec(vec A){
    vec B;
    B.initialize(A.rows,A.cols);
    for(int i=0;i<A.rows;i++){
        for(int j=0;j<A.cols;j++){
            if (A[i][j] < 0.0001){
                B[i][j]=0;
            }else{
                B[i][j]=A[i][j];
            }
        }
    }
    return B;
}
int baum_welch(vec A, vec B, vec p,list obsv,int maxiter, double acc){
    tuple <vec, vec, vec> out = estimation(A, B, p, obsv);
    vec newA = get<0>(out);
    vec newB = get<1>(out);
    vec newp = get<2>(out);
    vec tmp;
    int i=0;
    double error = 100;
    while(i<maxiter && error > acc){
        out = estimation(newA, newB, newp, obsv);
        tmp = get<0>(out);
        error = diff(newA,tmp);
        //cout<< "error :" << error << endl;;
        newA = fix_vec(get<0>(out));
        newB = fix_vec(get<1>(out));
        newp = fix_vec(get<2>(out));
        i++;
    }
    print_array_inline(newA);
    cout << endl;
    print_array_inline(newB);
    //print_vec(newA);
    
    
}

/////
// C
/////

///////////
// DUCKHUNT
///////////

/////////////
// Programme
/////////////

int main (){
    vec A = getinput();
    vec B = getinput();
    vec p = getinput();
    list obsv = get_observations();
    vec c;
    c.initialize(obsv.size,1);
    return 0;
    //print_vec(A);
    //print_vec(B);
    //print_vec(p);
    //print_list(obsv);
    /*
    pair <vec, vec> out = a_pass_norm(obsv,A,B,p);
    for(int i=0;i<5;i++){
        print_vec(extrct_row(out.first,i));
    }
    for(int i=0;i<5;i++){
        print_vec(extrct_row(out.second,i));
    }
    */
    /*
    print_vec(extrct_row(out.first,0));
    cout << endl;
    double sum_alpha = 0;
    for(int i=0; i<A.rows; i++){
        sum_alpha += Alpha.pi[obsv.size-1][i];
    }
    */
    //cout<< sum_alpha;
    /*
    /*
    cout << endl;
    printf("A");
    cout << endl;
    print_vec(A);

    cout << endl;
    printf("B");
    cout << endl;
    print_vec(B);
    
    cout << endl;
    printf("pi");
    cout << endl;
    print_vec(p);

    cout << endl;
    printf("Observations");
    cout << endl;
    print_list(obsv);
    cout << endl;
    */
    /*
    cout << endl;
    printf("output file");
    cout << endl;
    vec tmp = mat_mul(p,A);
    print_array_inline(mat_mul(tmp,B));
    
    cout << endl;
    printf("Observations");
    cout << endl;
    print_list(obsv);

    cout << endl;
    printf("Alpha t");
    cout << endl;
    print_vec(Alpha);
    */
    //pair<vec, vec> out = maximus(A);
    //print_vec(out.first);
    //print_vec(out.second);
    //cout << sum_alpha;
    //print_vec(Viterbi(A,B,p,obsv));
    //vec test;
    //multi_vec m_test (1,3,3);
    //m_test[0].initialize(3,3,3);
    //m_test[0][1][2] = 10;
    //cout << m.test[0]
    //vec X;
    //X.initialize(3,3,2.5);
    //vec Y;
    //Y.initialize(3,3,3);
    //X=X/2;
    //print_vec(X);
    //obsv.size = 1000; //selected observations !!!!!!!!!!
    //tuple <vec, vec, vec> out = estimation(A, B, p, obsv);
    //vec newA= get<0>(out);
    //vec newB = get<1>(out);
    //vec newpi = get<2>(out);
    //print_vec(newA);
    //print_vec(newB);
    //print_vec(newpi);
    //pair <vec, vec> outalpha= a_pass_norm(obsv,A,B,p);
    //vec Beta = b_pass_norm(A, B, obsv, obsv.size,outalpha.second);
    //print_vec(Beta);
    baum_welch(A, B, p, obsv, 50, 1e-04);
    return 0;
}
