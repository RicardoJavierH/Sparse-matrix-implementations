//
// Created by Ricardo on 24/06/2023.
//

#include "libreria.h"
#include "sparsemat.h"
#include "sparsevec.h"
#include <vector>
#include <iostream>

int main(){
    std::vector<std::vector<double>> M = {
           { 1., 0., 2., 0., 0. },
           { 0., 3., 0., 0., 0. },
           { 0., 0., 4., 0., 5. },
           { 6., 0., 0., 0., 0. },
           { 0., 7., 0., 8., 0. },
    };

    realVec B = {0., 0., 5., -1., 0., 0., 1.};
    realVec B2 = {0., 1., 5., 0., 0., 1., 0.};

    int nrows = 5;
    realVec B3 = {1., 1., 1., 1., 1.};
    realVec valA;
    intVec JA;
    intVec IA = {0};
    esparcifica(M,valA,IA,JA);
    printFullMatrix(M);
    printFullVector(valA, (char*)"valA = ");
    printFullVector(IA, (char*)"IA = ");
    printFullVector(JA, (char*)"JA = ");

    realVec out;
    //spMatrixVectorProd(nrows, valA, IA, JA, B3, out);
    //printVector(out, (char*)"AB = ");

    //******* Sparse vector using class data structure ********
    std::cout << "*** Sparse vector Implementation with class structure ***" << std::endl;

    VecSparse spB(B);
    spB.PrintJA((char*)"JB=");
    spB.PrintValA((char*)"ValB=");

    VecSparse spB2(B2);
    spB2.PrintJA((char*)"JB2=");
    spB2.PrintValA((char*)"ValB2=");

    VecSparse* spC;
    symbolicSpVecVecSum(spB,spB2,spC);
    NumericalSpVecVecSum(spB,spB2,spC);
    spC->PrintJA((char*)"JspC");
    spC->PrintValA((char*)"JspC");

    //******* Sparse matrix using class data structure ********
    std::cout << "*** Sparse matrix Implementation with class structure ***" << std::endl;
    MatSparseCSR sparseA(M);
    MatSparseCSR sparseB(M);
    sparseA.PrintValA((char*)"valA=");
    sparseA.PrintIA((char*)"IA=");
    sparseA.PrintJA((char*)"JA=");
    realVec resultVec(nrows);
    vectorSpMatProduct(B3, sparseA, resultVec);
    printFullVector(resultVec, (char*)"B*sparseA=\t");
    return 0;
}




