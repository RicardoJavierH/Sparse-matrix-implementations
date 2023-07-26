//
// Created Ricardo HP on 21/07/2023.
//
#include "sparsevec.h"

VecSparse::VecSparse(const realVec& fullVec) {
    size = fullVec.size();
    int i, j; //number of cols
    int nnz = 0;

    for (j = 0; j < size; j++) {
        if (fullVec[j] != 0) {
            valA.push_back(fullVec[j]); //agrega o inserta elementos al final del contenedor dinámico valA
            JA.push_back(j); //inserta j en el vector JA
            nnz++;
        }
    }

    Nnz = nnz;

}

VecSparse::VecSparse() {
    
}


void VecSparse::PrintValA(char* msg){
    std::cout << msg << "[ ";
    for_each(this->valA.begin(), this->valA.end(), [](int a) {
        std::cout << a << " ";
    });
    std::cout << "]" << std::endl;
}


void VecSparse::PrintJA(char* msg){
    std::cout << msg << "[ ";
    for_each(this->JA.begin(), this->JA.end(), [](int a) {
        std::cout << a << " ";
    });
    std::cout << "]" << std::endl;
}

int VecSparse::Size(){
    return this->size;
}

intVec* VecSparse::GetJA(){
    return &JA;
}

realVec* VecSparse::GetvalA(){
    return &valA;
}

void VecSparse::SetJA(const intVec& input){
    JA = input;
}
void VecSparse::SetValA(realVec& input){
    this->valA = input;
}
