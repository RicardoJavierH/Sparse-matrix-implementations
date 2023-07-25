//
// Created by Ricardo on 24/06/2023.
//

#ifndef LIBRERIA_H
#define LIBRERIA_H

#include <algorithm>
#include <vector>
#include "sparsemat.h"
#include "sparsevec.h"

typedef std::vector<double> realVec;
typedef std::vector<int> intVec;
typedef std::vector<std::vector<double>> realMat;
typedef std::vector<std::vector<int>> intMat;

template<class T>
void printFullMatrix(const std::vector<std::vector<T>>& M);

template<class T>
void printFullVector(const std::vector<T> & V, char* msg);

void esparcifica(const realMat& M, realVec& val, intVec& row_ptr, intVec& col_ind );

void spMatrixVectorProd(int n, const realVec& valA, const intVec& IA, const intVec& JA, const realVec& B, realVec& prod);

void vectorSpMatProduct(const realVec& vect, MatSparseCSR& mat, realVec& C );

void symbolicSpVecVecSum(VecSparse& A, VecSparse& B, VecSparse* out);
void NumericalSpVecVecSum(VecSparse& A, VecSparse& B, VecSparse* out);


MatSparseCSR spMatMatSymbolicSum(MatSparseCSR& A, MatSparseCSR& B);

#endif
