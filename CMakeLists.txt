cmake_minimum_required(VERSION 3.20)
project(parallel C)

FIND_PACKAGE(OpenMP REQUIRED)


if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -march=native -ffast-math ")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -march=native -ffast-math ")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS} ")
endif ()


set(CMAKE_CXX_STANDARD 17)
add_executable(main main.c)

TARGET_LINK_LIBRARIES(main -lpthread -lm)



