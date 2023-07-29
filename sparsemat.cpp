//
// Created by Ricardo on 24/06/2023.
//

#include "sparsemat.h"

MatSparseCSR::MatSparseCSR(const realMat& fullMat) {
    nrows = fullMat.size(); //number of rows
    ncols = fullMat[0].size();
    int i, j; //number of cols
    int nnz = 0;

    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ncols; j++) {
            if (fullMat[i][j] != 0) {
                valA.push_back(fullMat[i][j]); //agrega o incerta elementos al final del contenedor dinámico valA
                JA.push_back(j); //inserta j en el vector JA
                nnz++;
            }
        }
        IA.push_back(nnz);
    }
    Nnz = nnz;
}

MatSparseCSR::MatSparseCSR(){

}

void MatSparseCSR::PrintValA(char* msg){
    std::cout << msg << "[ ";
    for_each(this->valA.begin(), this->valA.end(), [](int a) {
        std::cout << a << " ";
    });
    std::cout << "]" << std::endl;
}

void MatSparseCSR::PrintIA(char* msg){
    std::cout << msg << "[ ";
    for_each(this->IA.begin(), this->IA.end(), [](int a) {
        std::cout << a << " ";
    });
    std::cout << "]" << std::endl;
}

void MatSparseCSR::PrintJA(char* msg){
    std::cout << msg << "[ ";
    for_each(this->JA.begin(), this->JA.end(), [](int a) {
        std::cout << a << " ";
    });
    std::cout << "]" << std::endl;
}

int MatSparseCSR::NCols(){
    return this->ncols;
}

int MatSparseCSR::NRows(){
    return this->nrows;
}

intVec* MatSparseCSR::GetIA(){
    return &IA;
}

intVec* MatSparseCSR::GetJA(){
    return &JA;
}

realVec* MatSparseCSR::GetValA(){
    return &valA;
}

void MatSparseCSR::SetIA(intVec& ia){
    this->IA = ia;
}
void MatSparseCSR::SetJA(intVec& ja){
    this->JA = ja;
}
void MatSparseCSR::SetValA(realVec vala){
    this->valA = vala;
}

void MatSparseCSR::SetNcols(int ncols){
    this->ncols = ncols;
}

void MatSparseCSR::SetNrows(int nrows){
    this->nrows = nrows;
}

void MatSparseCSR::SetNnz(int nnz){
    this->Nnz = nnz;
}
