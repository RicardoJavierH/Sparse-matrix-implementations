//
// Created by Ricardo on 24/06/2023.
//

#ifndef SPARSEMATRIX_SPARSEMAT_H
#define SPARSEMATRIX_SPARSEMAT_H

//#include "libreria.h"
#include <iostream>
#include <vector>
#include <algorithm>

typedef std::vector<double> realVec;
typedef std::vector<int> intVec;
typedef std::vector<std::vector<double>> realMat;
typedef std::vector<std::vector<int>> intMat;

class MatSparseCSR {
    public:
    MatSparseCSR(const realMat &fullMat); // constructor
    MatSparseCSR(); // constructor

    void PrintValA(char* msg);
    void PrintIA(char* msg);
    void PrintJA(char* msg);

    int NCols();
    int NRows();

    intVec* GetIA();
    intVec* GetJA();
    realVec* GetValA();

    void SetIA(intVec& ia);
    void SetJA(intVec& ja);
    void SetValA(realVec vala);
    void SetNcols(int ncols);
    void SetNrows(int nrows);
    void SetNnz(int nnz);

private:
    realVec valA;
    intVec IA={0};
    intVec JA;
    int ncols;
    int nrows = 0;
    int Nnz = 0; //number of non zeros
};


#endif //SPARSEMATRIX_SPARSEMAT_H
