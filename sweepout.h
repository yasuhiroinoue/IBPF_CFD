#ifndef ___SWEEPOUT_H
#define ___SWEEPOUT_H

template <class TEMPLATE> inline TEMPLATE myAbs(const TEMPLATE& a){
	if(a<0)return (-1.0*a);
	else return a;
}

template <class T> void sweep_out(const T* const in, T* const out, const int n){
    T* tmp_in;
    tmp_in = new T[n*n];
    memcpy(tmp_in, in, n*n*sizeof(T));    // Since we cannot change the matrix that was input.
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){            // Identity matrix.
            if(i==j) *(out+i*(n+1))=(T)1.0;
            else *(out+i+j*n)=(T)0.0;
        }
    }
    for(int i=0; i<n; i++){
        T p = *(tmp_in+i*(n+1));
        if(p==0){                        // When the head is 0.
            int j=0;
            for(j=i+1; j<n; j++){
                if(!(p = *(tmp_in+i+j*n)))break;
            }
            if(j==n) return;            // Rank deficiency. Should we return -1 here??
            else{                        // Row exchange.
                T* a;
                a = new T[n];
                for(int k=i; k<n; k++) a[k]=*(tmp_in+k+i*n);
                for(int k=i; k<n; k++) *(tmp_in+k+i*n)=*(tmp_in+k+j*n);
                for(int k=i; k<n; k++) *(tmp_in+k+j*n)=a[k];
                for(int k=0; k<n; k++) a[k]=*(out+k+i*n);
                for(int k=0; k<n; k++) *(out+k+i*n)=*(out+k+j*n);
                for(int k=0; k<n; k++) *(out+k+j*n)=a[k];
                delete[] a;
            }
        }
        if(p*p!=p){                    // Make the head 1.
            for(int j=i; j<n; j++) {*(tmp_in+j+i*n)/=p;}
            for(int j=0; j<n; j++) {*(out+j+i*n)/=p;}
        }
        for(int j=0; j<n; j++){        // Subtract from each column
            if(j==i) continue;
            if(p = *(tmp_in+i+j*n)){
                for(int k=i; k<n; k++) {*(tmp_in+k+j*n)-=(p**(tmp_in+k+i*n));}
                for(int k=0; k<n; k++) {*(out+k+j*n)-=(p**(out+k+i*n));}
            }
        }
    }
    delete[] tmp_in;
    return;
}

template <class T1,class T2> void sweep_out(T1 in[N_ARRAY][N_ARRAY],T2 out[N_ARRAY][N_ARRAY]){
    T2 tmp_in[N_ARRAY][N_ARRAY];
    for(int i=0;i<N_ARRAY;i++){
    for(int j=0;j<N_ARRAY;j++){                    // Because we cannot change the input matrix.
        tmp_in[i][j]=in[i][j];
    }}
    
    for(int i=0;i<N_ARRAY;i++){
    for(int j=0;j<N_ARRAY;j++){                                // Identity matrix.
        if(i==j)out[i][i]=(T2)1.0;
        else out[i][j]=0;
    }}
    for(int i=0;i<N_ARRAY;i++){
        T2 p = tmp_in[i][i];
        if(p==0){                                        // When the head is 0.
            int j=0;
            for(j=i+1;j<N_ARRAY;j++){
                if(!(p = tmp_in[j][i]))break;
            }
            if(j==N_ARRAY)return;                                // Rank fall. Should we return something like -1 here??
            else{                                        // Row exchange.
                T2 a[N_ARRAY];
                for(int k=i;k<N_ARRAY;k++)a[k]=tmp_in[i][k];
                for(int k=i;k<N_ARRAY;k++)tmp_in[i][k]=tmp_in[j][k];
                for(int k=i;k<N_ARRAY;k++)tmp_in[j][k]=a[k];
                for(int k=0;k<N_ARRAY;k++)a[k]=out[i][k];
                for(int k=0;k<N_ARRAY;k++)out[i][k]=out[j][k];
                for(int k=0;k<N_ARRAY;k++)out[j][k]=a[k];
            }
        }
        if(p*p!=p){                                        // Make the head 1.
            for(int j=i;j<N_ARRAY;j++){tmp_in[i][j]/=p;}
            for(int j=0;j<N_ARRAY;j++){out[i][j]/=p;}
        }
        for(int j=0;j<N_ARRAY;j++){                         // Subtract from each column.
            if(j==i)continue;
            if(p = tmp_in[j][i]){
                for(int k=i;k<N_ARRAY;k++){tmp_in[j][k]-=(p*tmp_in[i][k]);}
                for(int k=0;k<N_ARRAY;k++){out[j][k]-=(p*out[i][k]);}
            }
        }
    }
    return;
}
template <class T1,class T2,class T3> void lu(const T1 in[N_ARRAY][N_ARRAY],const T2 b[N_ARRAY], T3 x[N_ARRAY]){
    T2 tmp_in[N_ARRAY][N_ARRAY];
    for(int i=0;i<N_ARRAY;++i){
    for(int j=0;j<N_ARRAY;++j){                    // Because we cannot change the input matrix.
        tmp_in[i][j]=in[i][j];
    }}
    for(int i=0;i<N_ARRAY;++i){
        x[i]=b[i];
    }
    
    for(int i=0;i<N_ARRAY;++i){
        int d=i;
        double max = myAbs(tmp_in[i][i]);
        for(int j=i+1;j<N_ARRAY;++j){                    // Looking for the maximum absolute value.
            if(max<myAbs(tmp_in[j][i])){
                d=j;
                max=myAbs(tmp_in[j][i]);
            }
        }
        if(d!=i){                                        // Row exchange.
            T2 tmp=x[i];
            x[i]=x[d];
            x[d]=tmp;
            T2 a[N_ARRAY];
            for(int k=0;k<N_ARRAY;k++)a[k]=tmp_in[i][k];
            for(int k=0;k<N_ARRAY;k++)tmp_in[i][k]=tmp_in[d][k];
            for(int k=0;k<N_ARRAY;k++)tmp_in[d][k]=a[k];
        }
        if(tmp_in[i][i]){                                // Main content.
            for(int j=i+1;j<N_ARRAY;++j){
                tmp_in[j][i]/=tmp_in[i][i];
                for(int k=i+1;k<N_ARRAY;++k){
                    tmp_in[j][k]-=tmp_in[j][i]*tmp_in[i][k];
                }
            }
        }else {
            for(int j=0;j<N_ARRAY;++j)x[j]=0;
            return;                                     // Rank fall. Should we return something like -1 here??
        }
    }
    for(int i=0;i<N_ARRAY;++i){                         // Uy=b
        for(int j=0;j<i;++j){
            x[i]-=tmp_in[i][j]*x[j];
        }
    }
    for(int i=N_ARRAY-1;i>=0;--i){                      // Lx=y
        for(int j=N_ARRAY-1;j>i;--j){
            x[i]-=tmp_in[i][j]*x[j];
        }
        x[i]/=tmp_in[i][i];
    }
    
    return;
}
template <class T1,class T2> void lu_inv(const T1 in[N_ARRAY][N_ARRAY],T2 out[N_ARRAY][N_ARRAY]){
    T2 tmp_in[N_ARRAY][N_ARRAY];
    int index[N_ARRAY];
    for(int i=0;i<N_ARRAY;++i){
    for(int j=0;j<N_ARRAY;++j){                    // Because we cannot change the input matrix.
        tmp_in[i][j]=in[i][j];
    }}
    for(int i=0;i<N_ARRAY;++i){
        index[i]=i;
    }
    
    for(int i=0;i<N_ARRAY;++i){
        int d=i;
        double max = myAbs(tmp_in[i][i]);
        for(int j=i+1;j<N_ARRAY;++j){
            if(max<myAbs(tmp_in[j][i])){
                d=j;
                max=myAbs(tmp_in[j][i]);
            }
        }
        if(d!=i){
            T2 tmp=index[i];
            index[i]=index[d];
            index[d]=tmp;
            T2 a[N_ARRAY];
            for(int k=0;k<N_ARRAY;k++)a[k]=tmp_in[i][k];
            for(int k=0;k<N_ARRAY;k++)tmp_in[i][k]=tmp_in[d][k];
            for(int k=0;k<N_ARRAY;k++)tmp_in[d][k]=a[k];
        }
        if(tmp_in[i][i]){
            for(int j=i+1;j<N_ARRAY;++j){
                tmp_in[j][i]/=tmp_in[i][i];
                for(int k=i+1;k<N_ARRAY;++k){
                    tmp_in[j][k]-=tmp_in[j][i]*tmp_in[i][k];
                }
            }
        }else {
            return;
        }
    }
    T2 sum;
    for (int k=0;k<N_ARRAY;k++){
        for (int i=0;i<N_ARRAY;i++){
            if(index[i]==k){sum=1;}
            else {sum=0;}
            for(int j=0;j<i;j++){
                sum-=tmp_in[i][j]*out[j][k];
            }
            out[i][k]=sum;
        }
        for(int i=N_ARRAY-1;i>=0;i--){
            sum=out[i][k];
            for(int j=i+1;j<N_ARRAY;j++){
                sum-=tmp_in[i][j]*out[j][k];
            }
            out[i][k]=sum/tmp_in[i][i];
        }
    }
    return;
}
//#undef N_ARRAY





#endif

#if 0
//
// This comment is from a graduate named Mr. Deji.
// There's no particular need to include it, but I'm keeping it as a memento.
//
// This function calculates the inverse matrix by the sweep method.
// It assumes the array is an n*n matrix and loops accordingly,
// so it will probably throw a syntax error if used with a non-square matrix.
// It's safer to create a dedicated class.

// For the substitution method, see sweepout.cpp
// 1. Insert the head pointer of the array reduced to one dimension.
// 2. Push the head pointer of a static two-dimensional array.

// 3. If you have dynamically allocated a two-dimensional array with a double pointer and it is continuous,
//    insert the head pointer.
//    If you don't understand the meaning of array, array[0], &array[0], &array[0][0],
//    it's better not to try this.
// 4. If you have dynamically allocated a two-dimensional array with a double pointer and it is not continuous,
//    it's impossible.

// I'll study a bit more about specialization.

Usage:
#define N_ARRAY 100
#include "sweepout.h"
#undef N_ARRAY

Please include it like above.
N_ARRAY is the dimension of the matrix.

lu: When solving Ax=b for x,
    please substitute as lu(A,b,x).
sweep_out: When A_inv is the inverse matrix of A,
        please substitute as sweep_out(A,A_inv).
lu_inv:
        please substitute as lu_inv(A,A_inv).

If you are solving a system of linear equations directly,
you should compute in an area where a two-dimensional array can be statically allocated.
If more is needed, you should consider converging computations like SOR.
Please ask Ishida Sensei about SOR.

The result of a time attack is based on the calculation of a 100 by 100 matrix 1,000,000 times.
sweep_out 13141.9     3h39m01.9 s
lu         2237.67      37m17.67s
lu_inv    12652.8     3h30m52.8 s
Units are in seconds. Calculated on an old Opteron without compile options.
The computational order was supposed to be proportional to the cube of the degree. Just for reference.

Although I have directly calculated the determinant recursively,
I remember thinking it was not feasible to calculate it that way,
so I havent rewritten it.
I might need to create a function that solves equations without finding the inverse matrix by the sweep method,
but theoretically, the LU decomposition method is faster, and I dont have the energy to make it.


#endif

