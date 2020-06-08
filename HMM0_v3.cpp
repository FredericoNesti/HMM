hmm0_v3.cpp
Type
C++
Size
2 KB (1,904 bytes)
Storage used
0 bytes
Location
AI_A1
Owner
Theodor Panagiotakopoulos
Modified
Sep 5, 2019 by Theodor Panagiotakopoulos
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
void print_array_inline(array A){
    cout << A.rows << ' ' << A.cols;
    for (int i=0; i<A.rows; i++)    {
        for(int j=0; j< A.cols; j++){
            cout << ' ' << A.p[i][j];
        }
    }
    

}
array multiply_tables(array A, array B){
    array C;
    C.p = new float *[A.rows];
    C.rows = A.rows;
    C.cols = B.cols;
    for(int i=0;i<A.rows;i++){
        C.p[i] = new float[B.cols];
        for(int j=0;j<B.cols;j++){
            float sum = 0;
            for (int k=0;k<A.cols;k++){
                sum+= A.p[i][k]*B.p[k][j];
            }
            C.p[i][j] = sum;
        }
    }
    return C;
}
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
    float input[5000]={};
    getline(cin, line);
    //cout << line.length();
    int start = 0;
    int counter = 0;
    line += " ";
    bool flag = true;
    for (int i=0;i<line.length();i++){
        if (line.at(i) == ' ' && flag){
            input[counter] = stof(line.substr(start,i-start));
            start = i+1;
            counter++;
            flag = false;
        }else if (line.at(i) != ' '){
            flag = true;
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
    array tmp = multiply_tables(p,A);
    print_array_inline(multiply_tables(tmp,B));

return 0;
}
