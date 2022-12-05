cmake_minimum_required(VERSION 3.23)
project(Final_Project)

set(CMAKE_CXX_STANDARD 14)

if (CMAKE_SYSTEM_NAME MATCHES "Windows")
    include_directories(C:/TDM-GCC-64/include)

    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fopenmp")
    find_library(OPENMP_LIB libgomp.a C:\\TDM-GCC-64\\__installer\\downloaded\\lib\\lib\\gcc\\x86_64-w64-mingw32\\10.3.0)

elseif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    include_directories(/usr/local/include)

    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
    find_library(OPENMP_LIB libomp.dylib /usr/local/lib)
endif ()

add_executable(Final_Project
        main.cpp src/statistic_tool.cpp src/statistic_tool.h src/Payoff.h src/BSMModel.cpp src/BSMModel.h check_reduction.cpp check_asymptotic_result.cpp src/controlVariates.cpp src/controlVariates.h)

target_link_libraries(Final_Project LINK_PUBLIC ${OPENMP_LIB})
