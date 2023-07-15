//
// Created by Ricardo on 24/06/2023.
//

#include "sparsemat.h"
SparseCSR::SparseCSR(const Mat& fullMat) {
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

void SparseCSR::PrintValA(char* msg){
    std::cout << msg << "[ ";
    for_each(this->valA.begin(), this->valA.end(), [](int a) {
        std::cout << a << " ";
    });
    std::cout << "]" << std::endl;
}

void SparseCSR::PrintIA(char* msg){
    std::cout << msg << "[ ";
    for_each(this->IA.begin(), this->IA.end(), [](int a) {
        std::cout << a << " ";
    });
    std::cout << "]" << std::endl;
}

void SparseCSR::PrintJA(char* msg){
    std::cout << msg << "[ ";
    for_each(this->JA.begin(), this->JA.end(), [](int a) {
        std::cout << a << " ";
    });
    std::cout << "]" << std::endl;
}

int SparseCSR::NCols(){
    return this->ncols;
}

int SparseCSR::NRows(){
    return this->nrows;
}

Vec* SparseCSR::GetIA(){
    return &IA;
}

Vec* SparseCSR::GetJA(){
    return &JA;
}

Vec* SparseCSR::GetvalA(){
    return &valA;
}

void vecMatProduct(const Vec& vect, SparseCSR& mat, Vec& C ){
    int n = mat.NCols();
    int m = mat.NRows();
    Vec* ptria = mat.GetIA();
    Vec* ptrja = mat.GetJA();
    Vec* ptrvala = mat.GetJA();

     for(int irow=0; irow<m; irow++){
        for(int ii=0; ii<n; ii++){
            int IAA = (*ptria)[irow];
            int IAB = (*ptria)[irow+1]-1;
            
            if(IAB >= IAA){
                for(int k=IAA; k<IAB; k++){
                    int j = (*ptrja)[k];
                    C[j] = C[irow]+(*ptrvala)[k]*vect[irow];
                }
            }
        }
    }
}
