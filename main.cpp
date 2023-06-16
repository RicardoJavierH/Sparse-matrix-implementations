#include "libreria.h"

//Program to find CSR sparse matrix representation
//#include <algorithm>
//#include <iostream>
#include <vector>

//
int main()
{
   std::vector<std::vector<double>> M = {
           { 1., 0., 2., 0., 0. },
           { 0., 3., 0., 0., 0. },
           { 0., 0., 4., 0., 5. },
           { 6., 0., 0., 0., 0. },
           { 0., 7., 0., 8., 0. },
   };
    
   int nrows = 5;
   vec B = {1.,1.,1.,1.,1.};
   vec valA,JA;
   vec IA = {0};
   esparcifica(M,valA,IA,JA);
   printMatrix(M);
   printVector(valA, (char*)"valA = ");
   printVector(IA, (char*)"IA = ");
   printVector(JA, (char*)"JA = ");
   
    vec AB;
    matrixVectorProd(nrows, valA, IA, JA, B, AB);
    printVector(AB, (char*)"AB = ");
   return 0;
}




