cmake_minimum_required(VERSION 3.3)

project(cfe)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_COMPILER mpicxx)
# set(CMAKE_CXX_COMPILER mpiicpc)
# set(CMAKE_CXX_COMPILER g++)

include_directories(
    ~/include/fdps/src
    ~/workspace/src
    ${CMAKE_SOURCE_DIR}/include
)

add_compile_options(
    -ffast-math -funroll-loops
    -O3
    # -O0
    -g3
    # -fsanitize=address
    -DPARTICLE_SIMULATOR_THREAD_PARALLEL
    -DPARTICLE_SIMULATOR_MPI_PARALLEL
)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build/bin)

add_subdirectory(src)

add_executable(test.out ${CMAKE_SOURCE_DIR}/project/test/test.cc)
add_executable(speed.out ${CMAKE_SOURCE_DIR}/project/speed/speed.cc)
add_executable(debug.out ${CMAKE_SOURCE_DIR}/project/debug/debug.cc)
add_executable(plate.out ${CMAKE_SOURCE_DIR}/project/plate/plate.cc)
add_executable(relax.out ${CMAKE_SOURCE_DIR}/project/relax/relax.cc)
add_executable(linear.out ${CMAKE_SOURCE_DIR}/project/linear/linear.cc)
add_executable(tensile.out ${CMAKE_SOURCE_DIR}/project/tensile/tensile.cc)
add_executable(godunov.out ${CMAKE_SOURCE_DIR}/project/godunov/godunov.cc)
add_executable(torsion.out ${CMAKE_SOURCE_DIR}/project/torsion/torsion.cc)
add_executable(eigen_fn.out ${CMAKE_SOURCE_DIR}/project/eigen_fn/eigen_fn.cc)
add_executable(circ_linear.out ${CMAKE_SOURCE_DIR}/project/circ_linear/circ_linear.cc)

target_link_libraries(test.out PRIVATE ${MPI_LIBRARIES} _cfe)
target_link_libraries(speed.out PRIVATE ${MPI_LIBRARIES} _cfe)
target_link_libraries(debug.out PRIVATE ${MPI_LIBRARIES} _cfe)
target_link_libraries(plate.out PRIVATE ${MPI_LIBRARIES} _cfe)
target_link_libraries(relax.out PRIVATE ${MPI_LIBRARIES} _cfe)
target_link_libraries(linear.out PRIVATE ${MPI_LIBRARIES} _cfe)
target_link_libraries(tensile.out PRIVATE ${MPI_LIBRARIES} _cfe)
target_link_libraries(godunov.out PRIVATE ${MPI_LIBRARIES} _cfe)
target_link_libraries(torsion.out PRIVATE ${MPI_LIBRARIES} _cfe)
target_link_libraries(circ_linear.out PRIVATE ${MPI_LIBRARIES} _cfe)
