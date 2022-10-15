cmake_minimum_required(VERSION 3.23)
project(kmers)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(googletest)
include_directories(googletest/googletest/include)

add_executable(kmery
        #unittest.cpp
        main.cpp
        )

target_link_libraries(kmery gtest gtest_main)