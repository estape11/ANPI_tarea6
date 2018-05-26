cmake_minimum_required(VERSION 3.10)
project(ANPI_tarea6)

set(CMAKE_CXX_STANDARD 11)

enable_language(C)
enable_language(CXX)
enable_language(Fortran)
set(LAPACK_LIBRARIES lapack)
link_directories(${LAPACK_LIBDIRS})
add_executable(ANPI_tarea6 main.cpp)
target_link_libraries(ANPI_tarea6 ${LAPACK_LIBRARIES})
install(TARGETS ANPI_tarea6 DESTINATION ${CMAKE_INSTALL_PREFIX})
