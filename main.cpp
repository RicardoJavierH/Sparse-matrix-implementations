#include "libreria.h"

//Program to find sparse matrix representation using CSR
#include <algorithm>
#include <iostream>
#include <vector>

//
int main()
{
   std::vector<std::vector<int>> M = {
           { 1, 0, 2, 0, 0 },
           { 0, 3, 0, 0, 0 },
           { 0, 0, 4, 0, 5 },
           { 6, 0, 0, 0, 0 },
           { 0, 7, 0, 8, 0 },
   };

   vec B = {1,1,1,1,1};
   esparcifica(M); //papá

    
   return 0;
}




