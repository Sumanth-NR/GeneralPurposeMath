cmake_minimum_required(VERSION 3.21)
project(GeneralPurposeMath)

set(CMAKE_CXX_STANDARD 20)

add_subdirectory(lib)

add_executable(test test.cpp)
target_link_libraries(test sig)

add_executable(misc misc.cpp)
