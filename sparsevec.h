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

// Clase (plantilla) para construir objetos vector esparso
class VecSparse {
public:
    VecSparse(const realVec &fullVec); // constructor
    VecSparse(); // constructor

    void PrintValA(char* msg); // Método para imprimir valA
    void PrintJA(char* msg); // Método para imprimir JA

    int Size(); // Método para obtener el tamaño (del formato lleno) del vector esparso
    intVec* GetJA(); // Método para obtener el vector JA
    realVec* GetvalA(); // Método para obtener el vector valA

    void SetJA(const intVec& input); // Define JA
    void SetValA(realVec& input); // Define el valA

    void SetSize(int n);
    void SetNnz(int n);

        private:
    realVec valA;
    intVec JA;
    int size = 0;
    int Nnz = 0; //number of non zeros
};

#endif //SPARSEVEC_H
