//
// Created by Ricardo on 24/06/2023.
//

#include <iostream>
#include <vector>
#include "libreria.h"

//  Print full matrix
template <class T>
void printMatrix(const std::vector<std::vector<T>>& M)
{
   int m = M.size();
   int n = M[0].size();
   for (int i = 0; i < m; i++) {
       for (int j = 0; j < n; j++)
           std::cout << M[i][j] << " ";
       std::cout << std::endl;
   }
}


template<class T>
void printVector(const std::vector<T> & V, char* msg){
   std::cout << msg << "[ ";
   for_each(V.begin(), V.end(), [](int a) {
       std::cout << a << " ";
   });
   std::cout << "]" << std::endl;
}


// función que genera los tres vectores del formato CSR val, row_ptr, cold_ind para matrices cuadradas
void esparcifica(const realMat & M, realVec& val, intVec& row_ptr, intVec& col_ind ) //M vector de vectores
{
   int m = M.size(); //numero de filas
   int n = M[0].size(), i, j; //n = numero de columnas
//   std::vector<double> val;
//   std::vector<double> row_ptr = { 0 }; // row_ptr matríz de n+1 filas
//   std::vector<double> col_ind;
   int NNZ = 0;

   for (i = 0; i < m; i++) {
       for (j = 0; j < n; j++) {
           if (M[i][j] != 0) {
               val.push_back(M[i][j]); //agrega o incerta elementos al final de un contenedor dinámico en el vector val
               col_ind.push_back(j); //j inserta j en el vector col_ind

               // Count Number of Non Zeronumero de recuento distinto de cero
               // Elementos en la fila i
               NNZ++; // NNZ=CONT
           }
       }
       row_ptr.push_back(NNZ); //agrega elemento de nz al vector row_ptr
   }

//   printMatrix(M);
//   printVector(val, (char*)"val = ");
//   printVector(row_ptr, (char*)"row_ptr = ");
//   printVector(col_ind, (char*)"col_ind = ");
}

void matrixVectorProd(int n, const realVec& valA, const intVec& IA, const intVec& JA, const realVec& B, realVec& prod){
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

intVec symbolicSpVecVecSum(VecSparse& A, VecSparse& B){
    const int nA = A.Size();
    std::vector<int> IX(nA,0); // Crea IX con tamaño nA y entradas nulas

    const int nJA = A.GetJA()->size(); // Número de elementos de JA
    const int nJB = B.GetJA()->size(); // Número de elementos de JA
    std::vector<int> M;

    for(int i=0; i<nJA; i++){
        int j = (*A.GetJA())[i];
        IX[j] = 1;
        M.push_back(j); // Inserta j en el vector de unión M
    }
    for(int i=0; i<nJB; i++){
        int j = (*B.GetJA())[i];
        if(IX[j]==0){
            IX[j] = 1;
            M.push_back(j); // Inserta j en el vector de unión M
        }
    }
    return M;

}


MatSparseCSR spMatMatSymbolicSum(MatSparseCSR& A, MatSparseCSR& B){
    int N = A.NRows();
    int M = B.NCols();
    int IP = 1;
    intVec IX(M,0); // Create IX with size M and all entries zero
    intVec JC;
    intVec IC;

    for(int I=0; I<N; I++){
        IC[0]=1;
        int IAA = (*A.GetIA())[I];
        int IAB = (*A.GetIA())[I+1]-1;

        if(IAB >= IAA) {
            for (int JP = IAA; JP < IAB + 1; JP++) {
                int J = (*A.GetJA())[JP];
                JC[IP] = J;
                IP++;
                IX[J] = I;
            }

            int IBA = (*B.GetIA())[I];
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





