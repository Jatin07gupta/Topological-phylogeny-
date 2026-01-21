cmake_minimum_required(VERSION 3.10)
project(PhylogeneticTreeReconstruction LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# Find Eigen3
find_package(Eigen3 REQUIRED)

# Find OpenMP (optional - will work without it)
find_package(OpenMP)

# Find Eigen3
find_package(Eigen3 REQUIRED)

# Find OpenMP (optional - will work without it)
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found - parallel version will be enabled")
else()
    message(STATUS "OpenMP not found - using sequential version")
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${EIGEN3_INCLUDE_DIR})

# Source files
set(ALGORITHM2_SOURCES
    src/algorithm2.cpp
)

set(TEST_ALGORITHM2_SOURCES
    tests/test_algorithm2.cpp
    src/algorithm2.cpp
)

set(MAIN_ALGORITHM2_SOURCES
    src/main_algorithm2.cpp
    src/algorithm2.cpp
    src/data_loader.cpp
)

set(MLVMST_SOURCES
    src/mlvmst.cpp
    src/algorithm2.cpp
)

set(MAIN_MLVMST_SOURCES
    src/main_mlvmst.cpp
    src/mlvmst.cpp
    src/algorithm2.cpp
    src/data_loader.cpp
)

set(MLVMST_BORUVKA_SOURCES
    src/mlvmst_boruvka.cpp
    src/algorithm2.cpp
)

set(MAIN_MLVMST_BORUVKA_SOURCES
    src/main_mlvmst_boruvka.cpp
    src/mlvmst_boruvka.cpp
    src/mlvmst.cpp
    src/algorithm2.cpp
    src/data_loader.cpp
)

set(MLVMST_BORUVKA_NO_SORT_SOURCES
    src/mlvmst_boruvka_no_sort.cpp
    src/algorithm2.cpp
)

set(MAIN_MLVMST_BORUVKA_NO_SORT_SOURCES
    src/main_mlvmst_boruvka_no_sort.cpp
    src/mlvmst_boruvka_no_sort.cpp
    src/mlvmst_boruvka.cpp
    src/mlvmst.cpp
    src/algorithm2.cpp
    src/data_loader.cpp
)

set(MAIN_MLVMST_BORUVKA_PERFORMANCE_SOURCES
    src/main_mlvmst_boruvka_performance.cpp
    src/mlvmst_boruvka.cpp
    src/mlvmst_boruvka_no_sort.cpp
    src/mlvmst.cpp
    src/algorithm2.cpp
    src/data_loader.cpp
)

set(MLVMST_GPU_SOURCES
    src/mlvmst_gpu.cu
    src/algorithm2.cpp
)

set(MAIN_MLVMST_GPU_SOURCES
    src/main_mlvmst_gpu.cu
    src/mlvmst_gpu.cu
    src/mlvmst.cpp
    src/algorithm2.cpp
    src/data_loader.cpp
)

set(CLNJ_SOURCES
    src/clnj_fresh.cpp
)

set(MAIN_CLNJ_SOURCES
    src/main_clnj.cpp
    src/clnj_fresh.cpp
    src/mlvmst.cpp
    src/algorithm2.cpp
    src/data_loader.cpp
)

set(TEST_CLNJ_SOURCES
    test_clnj_direct.cpp
    src/clnj_fresh.cpp
)

set(MAIN_CLNJ_BORUVKA_SOURCES
    src/main_clnj_boruvka.cpp
    src/clnj_fresh.cpp
    src/mlvmst_boruvka.cpp
    src/mlvmst.cpp
    src/algorithm2.cpp
    src/data_loader.cpp
)

set(MAIN_GPU_MLVMST_TO_CLNJ_SOURCES
    src/main_gpu_mlvmst_to_clnj.cu
    src/mlvmst_gpu.cu
    src/clnj_fresh.cpp
    src/mlvmst.cpp
    src/algorithm2.cpp
    src/data_loader.cpp
)

set(DEBUG_ALGORITHM2_SOURCES
    debug_algorithm2.cpp
    src/algorithm2.cpp
    src/data_loader.cpp
)

# Build test executable
add_executable(test_algorithm2 ${TEST_ALGORITHM2_SOURCES})

# Build main executable
add_executable(algorithm2_main ${MAIN_ALGORITHM2_SOURCES})

# Build MLVMST executable
add_executable(mlvmst_main ${MAIN_MLVMST_SOURCES})

# Build MLVMST Borůvka executable
add_executable(mlvmst_boruvka_main ${MAIN_MLVMST_BORUVKA_SOURCES})

# Build MLVMST Borůvka No-Sort executable
add_executable(mlvmst_boruvka_no_sort_main ${MAIN_MLVMST_BORUVKA_NO_SORT_SOURCES})

# Build MLVMST Borůvka Performance comparison executable
add_executable(mlvmst_boruvka_performance_main ${MAIN_MLVMST_BORUVKA_PERFORMANCE_SOURCES})

# Build CLNJ executable
add_executable(clnj_main ${MAIN_CLNJ_SOURCES})

# Build CLNJ direct test executable
add_executable(test_clnj_direct ${TEST_CLNJ_SOURCES})

# Build CLNJ executable with Borůvka Parallel
add_executable(clnj_boruvka_main ${MAIN_CLNJ_BORUVKA_SOURCES})

# Build Algorithm 2 debug executable
add_executable(debug_algorithm2 ${DEBUG_ALGORITHM2_SOURCES})

# Build MLVMST GPU executable
add_executable(mlvmst_gpu_main ${MAIN_MLVMST_GPU_SOURCES})
set_target_properties(mlvmst_gpu_main PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)

# Build GPU MLVMST to CLNJ pipeline executable
add_executable(gpu_mlvmst_to_clnj_main ${MAIN_GPU_MLVMST_TO_CLNJ_SOURCES})
set_target_properties(gpu_mlvmst_to_clnj_main PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)

# Build flags
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-g -Wall -Wextra)
else()
    add_compile_options(-O3)
endif()

# Link OpenMP if available
if(OpenMP_CXX_FOUND)
    target_link_libraries(mlvmst_boruvka_main PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(mlvmst_boruvka_no_sort_main PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(mlvmst_boruvka_performance_main PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(clnj_boruvka_main PRIVATE OpenMP::OpenMP_CXX)
endif()

# Link CUDA libraries for GPU executables
find_library(CUDART_LIBRARY cudart)
if(CUDART_LIBRARY)
    target_link_libraries(mlvmst_gpu_main PRIVATE ${CUDART_LIBRARY})
    target_link_libraries(gpu_mlvmst_to_clnj_main PRIVATE ${CUDART_LIBRARY})
endif()

# Output directory
set_target_properties(test_algorithm2 PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(algorithm2_main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(mlvmst_main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(mlvmst_boruvka_main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(mlvmst_boruvka_no_sort_main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(mlvmst_boruvka_performance_main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(clnj_main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(clnj_boruvka_main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(test_clnj_direct PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(debug_algorithm2 PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(mlvmst_gpu_main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

set_target_properties(gpu_mlvmst_to_clnj_main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

