//
// Created by Ricardo on 24/06/2023.
//

#include "libreria.h"
#include "sparsemat.h"
#include "sparsevec.h"
#include <vector>
#include <iostream>

int main(){
    std::vector<std::vector<double>> M1 = {
           { 1., 0., 2., 0., 0. },
           { 0., 3., 0., 0., 0. },
           { 0., 0., 4., 0., 5. },
           { 6., 0., 0., 0., 0. },
           { 0., 7., 0., 8., 0. },
    };

    std::vector<std::vector<double>> M2 = {
            { 2., 0., 0., 1., 0. },
            { 0., 3., 0., 0., 0. },
            { 0., 0., 1., 1., 0. },
            { 0., 0., 1., 0., 2. },
            { 0., 1., 0., 0., 0. },
    };

    realVec v = {0., 0., 5., -1., 0., 0., 1.};
    realVec v2 = {0., 1., 5., 0., 0., 1., 0.};

    int nrows = 5;
    realVec B3 = {1., 1., 1., 1., 1.};
    realVec valA;
    intVec JA;
    intVec IA = {0};

    //******* Operaciones con vectores dispersos ********
    std::cout << "*** Sparse vector Implementation with class structure ***" << std::endl;

    VecSparse spV(v); // Crea el objeto spV (vector en formato esparso)
    spV.PrintJA((char*)"JB=");
    spV.PrintValA((char*)"ValB=");

    VecSparse spV2(v2);
    spV2.PrintJA((char*)"JB2=");
    spV2.PrintValA((char*)"ValB2=");

    VecSparse spSum;
    symbolicSpVecVecSum(spV,spV2,spSum);
    NumericalSpVecVecSum(spV,spV2,spSum);
    spSum.PrintJA((char*)"JspSum=");
    spSum.PrintValA((char*)"ValspSum=");

    //******* Operaciones con matrices dispersas ********
    std::cout << "*** Sparse matrix Implementation with class structure ***" << std::endl;
    MatSparseCSR spA(M1);
    spA.PrintValA((char*)"valspA=");
    spA.PrintIA((char*)"IspA=");
    spA.PrintJA((char*)"JspA=");
    int size = spA.NCols();
    realVec resultVec(size);
    vectorSpMatProduct(B3, spA, resultVec);
    printFullVector(resultVec, (char*)"B3*spA=");

    MatSparseCSR spB(M2);
    spB.PrintValA((char*)"valspB=");
    spB.PrintIA((char*)"IspB=");
    spB.PrintJA((char*)"JspB=");

    MatSparseCSR sum;
    spMatMatSymbolicSum(spA,spB,sum);
    sum.PrintIA((char*)"Isum=");
    sum.PrintJA((char*)"Jsum=");
    spMatMatNumericalSum(spA,spB,sum);
    sum.PrintValA((char*)"valsum=");
    
    MatSparseCSR prod;
    spMatMatSymbolicProd(spA,spB,prod);
    prod.PrintIA((char*)"Iprod=");
    prod.PrintJA((char*)"Jprod");
    spMatMatNumericalProd(spA, spB, prod);
    prod.PrintValA((char*)"valprod=");
    return 0;
}




