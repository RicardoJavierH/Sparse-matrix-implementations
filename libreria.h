#ifndef LIBRERIA_H
#define LIBRERIA_H

#include <vector>
typedef std::vector<double> vec;

void printMatrix(const std::vector<std::vector<double>> & M);

void printVector(const std::vector<double> & V, char* msg);

void esparcifica(const std::vector<std::vector<double>> & M, vec& val, vec& row_ptr, vec& col_ind );

void matrixVectorProd(int n, const vec& valA, const vec& IA, const vec& JA, const vec& B, vec& prod);








#endif
