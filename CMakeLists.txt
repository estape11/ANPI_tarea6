cmake_minimum_required(VERSION 3.10)
project(ANPI_tarea6)

set(CMAKE_CXX_STANDARD 11)

enable_language(C)
enable_language(CXX)
enable_language(Fortran)

find_package(LAPACK REQUIRED)

if(LAPACK_FOUND)
    message("Found Lapack")
    add_executable(ANPI_tarea6 main.cpp)
    target_link_libraries(ANPI_tarea6 lapack lapacke blas)
    message("Lapack link library is ${LAPACK_LIBRARIES}")
endif()