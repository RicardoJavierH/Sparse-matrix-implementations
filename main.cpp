
#include "libreria.h"
#include "sparsemat.h"
#include <vector>
#include <iostream>

int main()
{
    Mat M = {
           { 1., 0., 2., 0., 0. },
           { 0., 3., 0., 0., 0. },
           { 0., 0., 4., 0., 5. },
           { 6., 0., 0., 0., 0. },
           { 0., 7., 0., 8., 0. },
    };
    
    int nrows = 5;
    Vec B = {1., 1., 1., 1., 1.};
    Vec valA,JA;
    Vec IA = {0};
    esparcifica(M,valA,IA,JA);
    printMatrix(M);
    printVector(valA, (char*)"valA = ");
    printVector(IA, (char*)"IA = ");
    printVector(JA, (char*)"JA = ");

    Vec AB;
    matrixVectorProd(nrows, valA, IA, JA, B, AB);
    printVector(AB, (char*)"AB = ");

    //******* Sparse matrix using class data structure ********
    std::cout << "Implementacion con clases" << std::endl;
    SparseCSR sparseA(M);
    SparseCSR sparseB(M);
    sparseA.PrintValA((char*)"valA=");
    sparseA.PrintIA((char*)"IA=");
    sparseA.PrintJA((char*)"JA=");
    Vec prod(nrows);
    vecMatProduct(B, sparseA, prod);
    printVector(prod, (char*)"B*sparseA=\t");
    return 0;
}




