cmake_minimum_required(VERSION 3.10)
project(ANPI_tarea6)

set(CMAKE_CXX_STANDARD 11)

find_package(BLAS)
find_package(LAPACK)
if(LAPACK_FOUND AND BLAS_FOUND)
    set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

add_executable(ANPI_tarea6 main.cpp)
target_link_libraries(ANPI_tarea6 ${lapackblas_libraries})
