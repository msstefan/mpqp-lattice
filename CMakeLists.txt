cmake_minimum_required(VERSION 3.0.0)
project(mplattice VERSION 0.1.0)

# include(CTest)
# enable_testing()

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED)
# find_package (lcdd REQUIRED)
# find_package(OpenMP REQUIRED)

# add_subdirectory(include/eigen-cddlib)
add_executable(mplattice main.cpp src/eigen-cddlib/Polyhedron.cpp src/mpQP/mpQP.cpp src/Graph/Graph.cpp) # Also include src/Lattice/Lattice.cpp for lattice computation
target_include_directories(mplattice PUBLIC include/eigen-cddlib)
# target_include_directories(mplattice PUBLIC include/Lattice)
target_include_directories(mplattice PUBLIC include/mpQP)
target_include_directories(mplattice PUBLIC include/Graph)
target_include_directories(mplattice PUBLIC include/libigl/include)

target_link_libraries (mplattice Eigen3::Eigen)
target_link_libraries(mplattice cdd)
# target_link_libraries(mplattice OpenMP::OpenMP_CXX)
# target_link_libraries(mplattice eigen-cddlib)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(mplattice ${Boost_LIBRARIES})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION True) # Enable LTO


include(CPack)
