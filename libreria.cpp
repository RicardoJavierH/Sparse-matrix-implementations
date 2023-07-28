//
// Created by Ricardo on 24/06/2023.
//

#include <iostream>
#include <vector>
#include "libreria.h"

//  Print full matrix
template <class T>
void printFullMatrix(const std::vector<std::vector<T>>& M)
{
   int m = M.size();
   int n = M[0].size();
   for (int i = 0; i < m; i++) {
       for (int j = 0; j < n; j++)
           std::cout << M[i][j] << " ";
       std::cout << std::endl;
   }
}

template void printFullMatrix(const std::vector<std::vector<double>>& M);
template void printFullMatrix(const std::vector<std::vector<int>>& M);

template<class T>
void printFullVector(const std::vector<T> & V, char* msg){
   std::cout << msg << "[ ";
   for_each(V.begin(), V.end(), [](int a) {
       std::cout << a << " ";
   });
   std::cout << "]" << std::endl;
}

template void printFullVector(const std::vector<double> & V, char* msg);
template void printFullVector(const std::vector<int> & V, char* msg);

void spMatrixVectorProd(int n, const realVec& valA, const intVec& IA, const intVec& JA, const realVec& B, realVec& prod){
    for (int i=0; i< n; i++){
        double f = 0;
        double iaa = IA[i];
        double iab = IA[i+1]-1;
        
        if(iab >= iaa){
            for (int k=iaa; k<=iab; k++){
                f = f + valA[k]*B[JA[k]];
            }
        }
        prod.push_back(f);
        
    }
}

void vectorSpMatProduct(const realVec& vect, MatSparseCSR& mat, realVec& C ) {
    int n = mat.NCols();
    int m = mat.NRows();
    intVec *ptria = mat.GetIA();
    intVec *ptrja = mat.GetJA();
    realVec *ptrvala = mat.GetValA();

    for (int irow = 0; irow < m; irow++) {
        for (int ii = 0; ii < n; ii++) {
            int IAA = (*ptria)[irow];
            int IAB = (*ptria)[irow + 1] - 1;

            if (IAB >= IAA) {
                for (int k = IAA; k < IAB; k++) {
                    int j = (*ptrja)[k];
                    C[j] = C[irow] + (*ptrvala)[k] * vect[irow];
                }
            }
        }
    }
}

void symbolicSpVecVecSum(VecSparse& A, VecSparse& B, VecSparse& out){
    const int nA = A.Size();
    std::vector<int> IX(nA,0); // Crea IX con tamaño nA y entradas nulas

    const int nJA = A.GetJA()->size(); // Número de elementos de JA
    const int nJB = B.GetJA()->size(); // Número de elementos de JA
    std::vector<int> M;

    intVec JCaux;
    for(int i=0; i<nJA; i++){
        int j = (*A.GetJA())[i];
        IX[j] = 1;
        JCaux.push_back(j); // Inserta j en JCaux
    }
    for(int i=0; i<nJB; i++){
        int j = (*B.GetJA())[i];
        if(IX[j]==0){
            IX[j] = 1;
            JCaux.push_back(j); // Inserta j JCaux
        }
    }
    out.SetJA(JCaux);
    out.SetSize(nA);
    out.SetNnz(JCaux.size());
}

void NumericalSpVecVecSum(VecSparse& A, VecSparse& B, VecSparse& out) {
    const int nA = A.Size();
    realVec X(nA, 0); // Crea X con tamaño nA y entradas nulas

    const int nJA = A.GetJA()->size(); // Número de elementos de JA
    const int nJB = B.GetJA()->size(); // Número de elementos de JB

    for (int i = 0; i < nJA; i++) {
        int j = (*A.GetJA())[i];
        X[j] = (*A.GetvalA())[i];
    }

    for (int i = 0; i < nJB; i++) {
        int j = (*B.GetJA())[i];
        X[j] += (*B.GetvalA())[i];
    }

    intVec& JC = *out.GetJA();
    realVec ValC;
    for (int i=0; i < JC.size(); i++){
        int j = JC[i];
        ValC.push_back(X[j]);
    }
    out.SetValA(ValC);
}

void spMatMatSymbolicSum(MatSparseCSR& A, MatSparseCSR& B,MatSparseCSR& Sum){
    int N = A.NRows();
    int M = B.NCols();
    int IP = 1;
    intVec IX(M,0); // Create IX with size M and all entries zero
    intVec JC;
    intVec IC;

    for(int I=0; I<N; I++){
        IC.push_back(IP);
        int IAA = (*A.GetIA())[I];
        int IAB = (*A.GetIA())[I+1]-1;

        if(IAA <= IAB) {
            for (int JP = IAA; JP < IAB + 1; JP++) {
                int J = (*A.GetJA())[JP];
                JC[IP] = J;
                IP++;
                IX[J] = I;
            }

        }

        int IBA = (*B.GetIA())[I];
        int IBB = (*B.GetIA())[I+1]-1;

        if( IBA < IBB){
            for(int JP=IBA; JP<IBB+1; JP++){
                int J = (*B.GetJA())[JP];
                if(IX[J] != I){
                    JC[IP] = J;
                    IP++;
                }
            }
        }
    }
    IC[N+1] = IP;
}






