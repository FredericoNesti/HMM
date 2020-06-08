vec.h
Type
C
Size
5 KB (4,824 bytes)
Storage used
0 bytes
Location
version1
Owner
Theodor Panagiotakopoulos
Modified
Sep 15, 2019 by Theodor Panagiotakopoulos
Opened
11:43 PM by me
Created
Sep 15, 2019
Add a description
Viewers can download

#include <cstdlib>
#include <iostream>
#include<cmath>
class vec{
    public:
    int rows;
    int cols;
    double *data;
    //constructors
    vec(){
        rows = 0;
        cols = 0;
        data = 0;
    }
    vec(unsigned int n, unsigned int m){
        rows = n;
        cols = m;
        data = new double[n*m]();
    }
    vec(double * inarray){
        rows = (int) inarray[0];
        cols = (int) inarray[1];
        data = new double[rows*cols];
        for(int i=0;i<rows*cols;i++){
            data[i] = inarray[i+2];
        }
    }
    void recreate(unsigned int n, unsigned int m){
        rows = n;
        cols = m;
        if (data == 0){
            data = new double[n*m]();
        }else{
            delete[] data;
            data = new double[n*m]();
        }
        
    }
    vec(const vec& other){
        operator=(other);            
    }
    vec & operator=(const vec &other){
        rows = other.rows;
        cols = other.cols;
        data = new double[rows*cols];
        for(int i=0;i<rows*cols;i++){
            data[i] = other.data[i];
        }
        return *this;
    }
    ~vec(){
        //delete rows;
        //delete cols;
        if (data!=0){
            //std::cerr << "destructor called" << std::endl;
            delete[] data;
        }

    }
    //printing methods
    void print(){
        for(int i=0; i<rows; i++){
            for(int j=0; j<cols; j++){
                std::cerr << data[cols*i+j] << " ";
            }
            std::cerr << std::endl;
        }
    }
    void print_inline(){
        std::cout << rows << ' ' << cols;
        for (int i=0; i<rows; i++)    {
            for(int j=0; j<cols; j++){
                std::cout << ' ' << data[cols*i+j];
            }
        }
        std::cout << std::endl;
    }
    //get elements methods
    double & operator()(unsigned int i, unsigned int j){
        return data[cols*i+j];
    }
    vec get_row(int x){
        vec out(1,cols);
        for(int i=0;i<cols;i++){
            out.data[i] = data[cols*x+i];
        }
        return out;
    }
    vec get_col(int x){
        vec out(rows,1);
        for(int i=0;i<rows;i++){
            out.data[i] = data[cols*i +x];
        }
        return out;
    }
    // matrix multiplications

    vec operator*(const vec &B){ //matrix multiply
        vec C (rows, B.cols);
        double sum;
        for(int i=0;i<rows;i++){
            for(int j=0;j<B.cols;j++){
                sum=0;
                for(int ii=0;ii<cols;ii++){
                    sum+=data[cols*i+ii] * B.data[B.cols*ii + j];
                }
                C.data[B.cols*i+j] = sum;
            }
        }
        return C;
    }
    vec operator^(const vec &B){ //elementwise
        vec C (rows, cols);
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                C.data[cols*i+j] = data[cols*i+j]*B.data[cols*i+j];
            }
        }
        return C;
    }
    //math methods
    void randomize(){
        double x = (1.0/cols)/(cols-1)*0.9;
        for(int i=0;i<rows;i++){
            double sum=0;
            for(int j=0;j<cols-1;j++){
                double rnd = random()%10 -6;
                data[i*cols +j] = (1.0/cols)+rnd*x;
                sum += data[i*cols +j];
            }
            data[i*cols + cols-1] = 1-sum;
        }
    }
    vec trans(){
        vec C (cols,rows);
        for(int i=0;i<cols;i++){
            for(int j=0;j<rows;j++){
                C.data[i*rows + j] = data[j*cols + i];
            }
        }
        return C;
    }
    vec sum_cols(){
        vec c (rows,1);
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                c.data[i] += data[i*cols + j];
            }
        }
        return c;
    }
    std::pair<vec, vec> maximus(){
        vec out_pr (rows, 1);
        vec out_indx( rows, 1);
        int bin;
        for(int i=0;i<rows;i++){
            out_pr(i,0) = data[i*cols];
            bin= 1;
            out_indx(i,0) = 1;
            for(int j=1;j<cols;j++){
                bin*=2;
                if(out_pr(i,0) < data[i*cols+j]){
                    out_pr(i,0) = data[i*cols+j];
                    out_indx(i,0) = bin;
                }else if(out_pr(i,0) == data[i*cols+j]){
                    out_indx(i,0) += bin;
                }
            }
        }
        std::pair <vec, vec> OUT (out_pr,out_indx);
        return OUT;
    }
        // void methods
    void operator+=(const vec &B){
        for(int i=0;i<cols*rows;i++){
            data[i]+=B.data[i];
        }
    }
    void operator/=(const double x){
        for(int i=0;i<cols*rows;i++){
            data[i]/=x;
        }
    }
    void insert_row(const vec &A, int loc){
        for(int i=0;i<cols;i++){
            data[loc*cols + i] = A.data[i];
        }
    }
};
