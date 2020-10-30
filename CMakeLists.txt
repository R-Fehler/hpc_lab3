cmake_minimum_required(VERSION 3.10)
project(lab3 C)

set(CMAKE_C_STANDARD 11)

 add_executable(gameoflife_mpi gameoflife_mpi.c)
add_executable(helloworld helloworld.c)
add_executable(phase-field-mpi phase-field-mpi.c)

# if(MPI_FOUND)
# include_directories(SYSTEM ${MPI_INCLUDE_PATH})
# target_link_libraries(gameoflife_mpi ${MPI_C_LIBRARIES})


# endif()

find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(MPI_INCLUDE_PATH)
target_link_libraries(helloworld ${MPI_LIBRARIES})
target_link_libraries(gameoflife_mpi  ${MPI_LIBRARIES})
target_link_libraries(phase-field-mpi  ${MPI_LIBRARIES})