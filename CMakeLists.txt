cmake_minimum_required(VERSION 3.23)
project(kmers)

set(CMAKE_CXX_STANDARD 17)

add_subdirectory(googletest)
include_directories(googletest/googletest/include)

add_executable(kmers
        unittest.cpp
        #main.cpp
        )

target_link_libraries(kmers gtest gtest_main)
