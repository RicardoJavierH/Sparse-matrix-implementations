#include <algorithm>
#include <iostream>
#include <vector>
#include "libreria.h"

//  función para poder imprimir para matrix rectangular
void printMatrix(const std::vector<std::vector<int>> & M)
{
   int m = M.size();
   int n = M[0].size();
   for (int i = 0; i < m; i++) {
       for (int j = 0; j < n; j++)
           std::cout << M[i][j] << " ";
       std::cout << std::endl;
   }
}


// para imprimir los vectores val, row_ptr, col_ind.
void printVector(const std::vector<int> & V, char* msg)
//define una función llamada printVector que toma un vector de enteros y un mensaje como parámetros,
// e imprime el vector junto con el mensaje, sin realizar ninguna modificación en el vector original.
{


   std::cout << msg << "[ ";
   for_each(V.begin(), V.end(), [](int a) {
       std::cout << a << " ";
   });
   std::cout << "]" << std::endl;
}


// función que genera los tres vectores del formato CSR val, row_ptr, cold_ind para matrices cuadradas
void esparcifica(const std::vector<std::vector<int>> & M) //M vector de vectores
{
   int m = M.size(); //numero de filas
   int n = M[0].size(), i, j; //n = numero de columnas
   std::vector<int> val;
   std::vector<int> row_ptr = { 0 }; // row_ptr matrís de n+1 filas
   std::vector<int> col_ind;
   int NNZ = 0;


   for (i = 0; i < m; i++) {
       for (j = 0; j < n; j++) {
           if (M[i][j] != 0) {
               val.push_back(M[i][j]); //agrega o incerta elementos al final de un contenedor dinámico en el vector val
               col_ind.push_back(j); //j inserta j en el vector col_ind


               // Count Number of Non Zeronumero de recuento distinto de cero
               // Elementos en la fila i
               NNZ++; // NNZ=CONT
           }
       }
       row_ptr.push_back(NNZ); //agrega elemento de nz al vector row_ptr
   }


   printMatrix(M);
   printVector(val, (char*)"val = ");
   printVector(row_ptr, (char*)"row_ptr = ");
   printVector(col_ind, (char*)"col_ind = ");
}

void matrixVectorProd(int n, const vec& valA, const vec& IA, const vec& JA, const vec& B, vec& C){
    for (int i=0; i< n; i++){
        int f = 0;
        double iaa = IA[i];
        double iab = IA[i+1]-1;
        
        if(iab >= iaa){
            for (int k=iaa; k<=iab; k++){
                f = f + valA[k]*B[JA[k]];
            }
        }
        C[i] = f;
        
    }
}

