//
// Created by Ricardo on 24/06/2023.
#ifndef SPARSEMATRIX_SPARSEMAT_H
#define SPARSEMATRIX_SPARSEMAT_H

#include "libreria.h"
#include <iostream>
//#include <vector>
class sparseCSR {
    public:
    sparseCSR(const mat &fullMat);
    void printValA(char* msg);
    void printIA(char* msg);
    void printJA(char* msg);

    private:


    vec valA;
    vec IA;
    vec JA;
};

#endif //SPARSEMATRIX_SPARSEMAT_H
