hmm0_v2.cpp
Type
C++
Size
2 KB (2,363 bytes)
Storage used
0 bytes
Location
AI_A1
Owner
Theodor Panagiotakopoulos
Modified
Sep 5, 2019 by Theodor Panagiotakopoulos
Opened
11:40 PM by me
Created
Sep 5, 2019
Add a description
Viewers can download

#include <cstdlib>
#include <iostream>
#include <string.h>
using namespace std;
struct array{
    int rows;
    int cols;
    float **p;
};
/*
float **multiply_tables(float **A, float **B){
    float ** C = new float *[A.size()];
        for(int i=0;i<A.size();i++){
            C[i] = new float[B[0].size()];
            for(int j=0;j<B[0].size();j++){
                double sum = 0;
                for (int k=0;k<A[i].size();k++){
                    sum+= A[i][k]*B[k][j];
                }
                C[i][j] = sum;
            }
        }
        return C;
}

int array_size(float A[]){
    return A.size();
}
*/
void print_array(array A){
    for(int i=0;i<A.rows;i++){
        for(int j=0;j<A.cols;j++){
            cout << A.p[i][j] << ' ';
        }
        cout << '\n';
    }
}
array getinput(){
    string line;
    float input[300]={};
    getline(cin, line);
    //cout << line.length();
    int start = 0;
    int counter = 0;
    line += " ";
    for (int i=0;i<line.length();i++){
        if (line.at(i) == ' '){
            input[counter] = stof(line.substr(start,i-start));
            start = i+1;
            counter++;
        }
    }
    array A;
    A.rows = (int)input[0];
    A.cols = (int)input[1];
    A.p = new float *[A.rows];
    for(int i=0;i<A.rows;i++){
        A.p[i] = new float[A.cols];
        for(int j=0;j<A.cols;j++){
            A.p[i][j] = input[2 + i*A.cols+ j];
        }
    }
    return A;
}

int main (){
    array A = getinput();
    array B = getinput();
    array p = getinput();
    //float **B = getinput();
    //float **p = getinput();
    //float **tmp = multiply_tables(p,A);
    //print_table(A);
    //cout << A.rows;
    print_array(A);
    print_array(B);
    print_array(p);



    //printf("\n%s",line.c_str());
    //cout << input[3];
    return 0;
   
    //float input[1000] = {};
    /*float n;
    int i=0;
    while((scanf("%f",&n)) == 1)
    {
        printf("%f",n);
        i++;
    }
    
    return 0;*/
    //cin >> input;
    //cout << sizeof(input);
    /*for (int i=0;i<strlen(input);i++){
        cout << input[i] << "\n";
    }
    return 0;
    
    float A[input[0]][input[2]];
    float x[100];
    for (int i=0;i<input[0];i++){
        for (int j=0;j<input[2];j++){
            x[input[0]*i+j] = input[4+2*(input[0]*i + j)]; 
        }
    }
    cout << &x;*/
return 0;
}
