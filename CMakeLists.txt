cmake_minimum_required(VERSION 3.10)
project(Lab_1)

set(CMAKE_CXX_STANDARD 17)

add_executable(Lab_1 project/src/main.cpp project/src/Matrix.cpp
        project/include/Matrix.h project/src/LinearSystems.cpp project/include/LinearSystems.h
        project/src/Interpolation.cpp project/include/Interpolation.h project/src/EigenValues.cpp
        project/include/EigenValues.h project/src/NonlinearSolve.cpp project/include/NonlinearSolve.h)