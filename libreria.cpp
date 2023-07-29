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
    int ncol = mat.NCols();
    int nrow = mat.NRows();
    int size = vect.size();
    if (size != nrow){
        std::cerr << "The dimensions aren't compatible!!" << std::endl;
        exit(1);
    }
    
    intVec *ptria = mat.GetIA();
    intVec *ptrja = mat.GetJA();
    realVec *ptrvala = mat.GetValA();

    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            int IAA = (*ptria)[i];// Inicio de la fila i en JA y valA
            int IAB = (*ptria)[i + 1] - 1; // fin de la fila i en JA y valA

            if (IAB >= IAA) {
                for (int k = IAA; k < IAB; k++) {
                    int J = (*ptrja)[k];
                    C[J] = C[i] + (*ptrvala)[k] * vect[i];
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
    int nrows = A.NRows();
    int ncols = B.NCols();
    intVec IX(ncols,0); // Create IX with size M and all entries zero
    intVec JC;
    intVec IC={0};
    int IP=0;
    
    for(int I=0; I<nrows; I++){
        int IAA = (*A.GetIA())[I]; //Inicio de la fila I en JA y valA
        int IAB = (*A.GetIA())[I+1]-1; //final de la fila I en JA y valA

        //Inserta en JC los índices de las columnas no nulas de la fila I en A
        // y en XI son marcadas con I+1 las posiciones no nulas
        if(IAA <= IAB) {
            for (int JP = IAA; JP < IAB + 1; JP++) {
                int J = (*A.GetJA())[JP];
                JC.push_back(J);
                IX[J]=I+1;
                IP++;
            }

        }

        int IBA = (*B.GetIA())[I];//Inicio de la fila I en JB y valB
        int IBB = (*B.GetIA())[I+1]-1;//Fin de la fila I en JB y valB

        //Inserta en JC los índices de las columnas no nulas de la fila I en B
        if( IBA < IBB){
            for(int JP=IBA; JP<IBB+1; JP++){
                int J = (*B.GetJA())[JP];
                if(IX[J] != I+1){
                    JC.push_back(J);
                    IP++;
                }
            }
        }
        
        IC.push_back(IP);
        
    }
    
    Sum.SetIA(IC);
    Sum.SetJA(JC);
}

void spMatMatNumericalSum(MatSparseCSR& A, MatSparseCSR& B, MatSparseCSR& Sum){
    int nrows = A.NRows();
    int ncols = B.NCols();
    realVec X(ncols,0); // Create IX with size M and all entries zero
    intVec& IC = *Sum.GetIA();
    intVec& JC = *Sum.GetJA();
    realVec valC(JC.size());
    intVec& IA = *A.GetIA();
    intVec& JA = *A.GetJA();
    realVec& valA = *A.GetValA();
    intVec& IB = *B.GetIA();
    intVec& JB = *B.GetJA();
    realVec& valB = *B.GetValA();

    for(int I=0; I<nrows; I++){
        int IH = I+1;
        int ICA = IC[I];
        int ICB = IC[I+1]-1;
        
        if (ICA <= ICB){
            for (int IP=ICA; IP<ICB+1; IP++){
                X[JC[IP]] = 0;
            }
            
            int IAA = IA[I];
            int IAB = IA[I+1]-1;
            
            if (IAA <= IAB){
                for (int IP=IAA; IP<IAB+1; IP++){
                    X[JA[IP]] = valA[IP];
                }
            }
            
            int IBA = IB[I];
            int IBB = IB[I+1]-1;
            
            if (IBA <= IBB){
                for (int IP=IBA; IP<IBB+1; IP++){
                    int J = JB[IP];
                    X[J] = X[J] + valB[IP];
                }
                
                for (int IP=ICA; IP<ICB+1;IP++){
                    valC[IP] = X[JC[IP]];
                }
                
            }
            
        }
        
    }
    
    Sum.SetValA(valC);
    
}

void spMatMatSymbolicProd(MatSparseCSR& A, MatSparseCSR& B, MatSparseCSR& prod){
    if (A.NCols() != B.NRows()){
        std::cerr << "Incorrect dimensions for Matrix product!! "<<std::endl;
        exit(1);
    }
    
    intVec& IA = *A.GetIA();
    intVec& JA = *A.GetJA();
    intVec& IB = *B.GetIA();
    intVec& JB = *B.GetJA();
    int NP = A.NRows();
    int NQ = A.NCols();
    int NR = B.NCols();

    prod.SetNrows(NP);
    prod.SetNcols(NR);
    
    intVec IC={0};
    intVec JC;
    intVec IX(NR); // Multiple switch

    int IP = 0;
    for (int I=0; I<NP; I++){
        int IAA = IA[I];//Inicio de la fila I en JA y val A
        int IAB = IA[I+1]-1;//Final de la fila I en JA y val A
        
        if (IAA <= IAB){
            for (int JP=IAA; JP<IAB+1; JP++){
                int J = JA[JP];
                int IBA = IB[J];//Inicio de la fila I en JA y val A
                int IBB = IB[J+1]-1;//Final de la fila I en JA y val A
                if (IBA <= IBB){
                    for (int KP=IBA; KP<IBB+1; KP++){
                        int K = JB[KP];
                        if (IX[K] != I+1){
                            JC.push_back(K);
                            IP++;
                            IX[K]=I+1;
                        }
                    }
                }

            }
        }
        
        IC.push_back(IP);
        
    }
    IC[NP+1]=IP;
    prod.SetIA(IC);
    prod.SetJA(JC);
}

void spMatMatNumericalProd(MatSparseCSR& A, MatSparseCSR& B, MatSparseCSR& prod){
    intVec& IA = *A.GetIA();
    intVec& JA = *A.GetJA();
    realVec& valA = *A.GetValA();
    intVec& IB = *B.GetIA();
    intVec& JB = *B.GetJA();
    realVec& valB = *B.GetValA();
    intVec& IC = *prod.GetIA();
    intVec& JC = *prod.GetJA();
    realVec& valC = *prod.GetValA();
    
    int NP = A.NRows();
    realVec X(prod.NCols());
    
    for (int I=0; I<NP; I++){
        int ICA = IC[I];// Inicio de la fila I en JC y val C
        int ICB = IC[I+1]-1;//Final de la fila I en JC y valC
        
        if (ICA <= ICB){
            for (int J = ICA; J<ICB+1; J++){
                X[JC[J]] = 0.;
            }
            
            int IAA = IA[I];
            int IAB = IA[I+1]-1;
            
            for (int JP=IAA; JP<IAB+1; JP++){
                int J = JA[JP];
                double A = valA[JP];
                
                int IBA = IB[J];
                int IBB = IB[J+1]-1;
                
                if (IBA <= IBB){
                    for (int KP=IBA; KP<IBB+1; KP++){
                        int K = JB[KP];
                        X[K] = X[K]+A*valB[KP];
                    }
                }
            }
        
            for (int J=ICA; J<ICB+1; J++){
                valC.push_back(X[JC[J]]);
            }
            
        }
    }
    prod.SetValA(valC);
}


