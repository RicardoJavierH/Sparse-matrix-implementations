//
// Created by Ricardo on 24/06/2023.
//

#include "sparsemat.h"
sparseCSR::sparseCSR(const mat& fullMat) {
    int m = fullMat.size(); //número de filas
    int n = fullMat[0].size(), i, j; //número de columnas
    int NNZ = 0;

    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if (fullMat[i][j] != 0) {
                valA.push_back(fullMat[i][j]); //agrega o incerta elementos al final de un contenedor dinámico en el vector val
                JA.push_back(j); //j inserta j en el vector col_ind

                // Count Number of Non Zeronumero de recuento distinto de cero
                // Elementos en la fila i
                NNZ++; // NNZ=CONT
            }
        }
        IA.push_back(NNZ); //agrega elemento de nz al vector row_ptr
    }
}

void sparseCSR::printValA(char* msg){
    std::cout << msg << "[ ";
    for_each(this->valA.begin(), this->valA.end(), [](int a) {
        std::cout << a << " ";
    });
    std::cout << "]" << std::endl;
}

void sparseCSR::printIA(char* msg){
    std::cout << msg << "[ ";
    for_each(this->IA.begin(), this->IA.end(), [](int a) {
        std::cout << a << " ";
    });
    std::cout << "]" << std::endl;
}

void sparseCSR::printJA(char* msg){
    std::cout << msg << "[ ";
    for_each(this->JA.begin(), this->JA.end(), [](int a) {
        std::cout << a << " ";
    });
    std::cout << "]" << std::endl;
}