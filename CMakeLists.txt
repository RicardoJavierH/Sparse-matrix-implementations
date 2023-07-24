cmake_minimum_required(VERSION 3.12)
project(SparseMatrix)

set(CMAKE_CXX_STANDARD 17)

# Agrega el archivo.cpp adicional al proyecto
set(SOURCES
    libreria.cpp
    sparsemat.cpp
    sparsevec.cpp
)

# Agrega el archivo.h al proyecto
set(HEADERS
    libreria.h
    sparsemat.h 
    sparsevec.h
)

# Agrega el archivo main.cpp al proyecto
add_executable(SparseMatrix main.cpp ${SOURCES} ${HEADERS})
