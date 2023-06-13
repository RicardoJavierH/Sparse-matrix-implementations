#ifndef LIBRERIA_H
#define LIBRERIA_H

#include <vector>
typedef std::vector<int> vec;

void printMatrix(const std::vector<std::vector<int>> & M);

void printVector(const std::vector<int> & V, char* msg);

void esparcifica(const std::vector<std::vector<int>> & M);

void matrixVectorProd(int n, const vec& valA, const vec& IA, const vec& JA, const vec& B, vec& C);








#endif
