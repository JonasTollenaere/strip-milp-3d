project(strip-milp)

cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

### CPP ###
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
endif()

### MeshCore ###
add_compile_definitions(GLM_ENABLE_EXPERIMENTAL)
set(MESHCORE_BUILD_DEMOS OFF)
include(FetchContent)
FetchContent_Declare(MeshCore
        GIT_REPOSITORY git@github.com:JonasTollenaere/MeshCore.git
        GIT_TAG        v0.1.1
)
FetchContent_MakeAvailable(MeshCore)

file(GLOB SOURCES src/*.cpp)
include_directories(include)
add_compile_definitions(SOLUTION_DIR="${CMAKE_CURRENT_SOURCE_DIR}/solutions/")
add_compile_definitions(MODEL_DIR="${CMAKE_CURRENT_SOURCE_DIR}/models/")

#find_package(GUROBI REQUIRED) # FindGUROBI.cmake not found when using vcpkg toolchain, so we pass paths as explicit CMAKE variables
# Examples:
# set(GUROBI_INCLUDE_DIRS ${GUROBI_DIR}/include)
# set(GUROBI_CXX_LIBRARY ${GUROBI_DIR}/lib/gurobi_c++md2017.lib)
# set(GUROBI_CXX_DEBUG_LIBRARY ${GUROBI_DIR}/lib/gurobi_c++mdd2017.lib)
# set(GUROBI_LIBRARY ${GUROBI_DIR}/lib/gurobi110.lib)
if (DEFINED GUROBI_INCLUDE_DIRS AND DEFINED GUROBI_CXX_LIBRARY AND DEFINED GUROBI_CXX_DEBUG_LIBRARY AND DEFINED GUROBI_LIBRARY)

    add_executable(RunMILP targets/RunMILP.cpp ${SOURCES})
    target_include_directories(RunMILP PRIVATE ${GUROBI_INCLUDE_DIRS})
    target_link_libraries(RunMILP PRIVATE MeshCore)
    target_link_libraries(RunMILP PRIVATE optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})
    target_link_libraries(RunMILP PRIVATE ${GUROBI_LIBRARY})

    # The benchmark target does not require any UI dependencies (only MeshCoreBase) and could be build without QT
    add_executable(BenchmarkMILP targets/BenchmarkMILP.cpp ${SOURCES})
    target_include_directories(BenchmarkMILP PRIVATE ${GUROBI_INCLUDE_DIRS})
    target_link_libraries(BenchmarkMILP PRIVATE MeshCoreBase)
    target_link_libraries(BenchmarkMILP PRIVATE optimized ${GUROBI_CXX_LIBRARY} debug ${GUROBI_CXX_DEBUG_LIBRARY})
    target_link_libraries(BenchmarkMILP PRIVATE ${GUROBI_LIBRARY})

else ()
    message("GUROBI not found. Skipping Strip_Packing_MILP.exe")
endif ()

### VisualizeSlices ###
add_executable(VisualizeSlices targets/VisualizeSlices.cpp ${SOURCES})
target_link_libraries(VisualizeSlices PRIVATE MeshCore)

### VisualizeSolution ###
add_executable(VisualizeSolution targets/VisualizeSolution.cpp ${SOURCES})
target_link_libraries(VisualizeSolution PRIVATE MeshCore)