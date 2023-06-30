#ifndef LIBRERIA_H
#define LIBRERIA_H

#include <algorithm>
#include <vector>
typedef std::vector<double> Vec;
typedef std::vector<std::vector<double>> Mat;
void printMatrix(const Mat& M);

void printVector(const Vec & V, char* msg);

void esparcifica(const Mat& M, Vec& val, Vec& row_ptr, Vec& col_ind );

void matrixVectorProd(int n, const Vec& valA, const Vec& IA, const Vec& JA, const Vec& B, Vec& prod);




#endif
