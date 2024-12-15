cmake_minimum_required(VERSION 3.11 FATAL_ERROR)
project(BLTC LANGUAGES CXX CUDA)

add_executable(test BLTC.cu directsum.cu main.cpp)

set(CMAKE_CUDA_FLAGS "-arch=native")
