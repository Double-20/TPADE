cmake_minimum_required(VERSION 3.16.9)
project(TPADE)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(TPADE_CEC17
        main_CEC17.cpp
        Multimethod.cpp
        Multimethod.h
        Population.h
        CEC17Problems.cpp
        CEC17Problems.h)

add_executable(TPADE_CEC14
        main_CEC14.cpp
        Multimethod.cpp
        Multimethod.h
        Population.h
        CEC14Problems.cpp
        CEC14Problems.h)

add_executable(TPADE_CEC22
        main_CEC22.cpp
        Multimethod.cpp
        Multimethod.h
        Population.h
        CEC22Problems.cpp
        CEC22Problems.h)

