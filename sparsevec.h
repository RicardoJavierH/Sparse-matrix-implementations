//
// Created by Ricardo on 21/07/2023.
//

#ifndef SPARSEVEC_H
#define SPARSEVEC_H

//#include "libreria.h"
#include <iostream>
#include <vector>
#include <algorithm>

typedef std::vector<double> realVec;
typedef std::vector<int> intVec;

class VecSparse {
public:
    VecSparse(const realVec &fullVec);
    void PrintValA(char* msg);
    void PrintJA(char* msg);
    int Size();
    intVec* GetJA();
    realVec* GetvalA();
    void SetJA(intVec& input);
    void SetValA(realVec& input);

private:
    realVec valA;
    intVec JA;
    int size;
    int Nnz = 0; //number of non zeros
};

#endif //SPARSEVEC_H
