cmake_minimum_required(VERSION 3.23)
project(Final_Project)

set(CMAKE_CXX_STANDARD 14)

include_directories(/usr/local/include)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
find_library(OPENMP_LIB libomp.dylib /usr/local/lib)

add_executable(Final_Project
        main.cpp statistic_tool.cpp statistic_tool.h Payoff.h BSMModel.cpp BSMModel.h check_reduction.cpp check_asymptotic_result.cpp varianceReduction.cpp varianceReduction.h)

target_link_libraries(Final_Project LINK_PUBLIC ${OPENMP_LIB})