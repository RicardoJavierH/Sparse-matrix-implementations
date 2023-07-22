//
// Created by Ricardo on 21/07/2023.
//

#ifndef SPARSEVEC_H
#define SPARSEVEC_H

#include "libreria.h"
#include <iostream>

class VecSparse {
public:
    VecSparse(const Vec &fullVec);
    void PrintValA(char* msg);
    void PrintJA(char* msg);
    int Size();
    Vec* GetJA();
    Vec* GetvalA();

private:
    Vec valA;
    Vec JA;
    int size;
    int Nnz = 0; //number of non zeros
};

#endif //SPARSEVEC_H
