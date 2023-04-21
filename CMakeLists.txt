cmake_minimum_required(VERSION 3.20)

project(BS-SOLCTRA LANGUAGES C CXX)
option(VECTOR "Build for vector processor" OFF)
option(TREAT_WARNINGS_AS_ERRORS "Treat all warnings as errors" OFF)


set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -Wall -Wpedantic -Wextra")

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")

add_subdirectory(main)