#ifndef LIBRERIA_H
#define LIBRERIA_H

#include <algorithm>
#include <vector>
typedef std::vector<double> vec;
typedef std::vector<std::vector<double>> mat;
void printMatrix(const mat& M);

void printVector(const vec & V, char* msg);

void esparcifica(const mat& M, vec& val, vec& row_ptr, vec& col_ind );

void matrixVectorProd(int n, const vec& valA, const vec& IA, const vec& JA, const vec& B, vec& prod);








#endif
