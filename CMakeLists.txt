cmake_minimum_required(VERSION 3.12)
project(SparseMatrix)

set(CMAKE_CXX_STANDARD 11)

# Agrega el archivo.cpp adicional a tu proyecto
set(SOURCES
    libreria.cpp
)

# Agrega el archivo.h a tu proyecto
set(HEADERS
    libreria.h
)

# Agrega el archivo main.cpp a tu proyecto
add_executable(SparseMatrix main.cpp ${SOURCES} ${HEADERS})
