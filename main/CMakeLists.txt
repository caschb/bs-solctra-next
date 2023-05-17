set(EXECUTABLE_TARGET bs_solctra)

add_executable(${EXECUTABLE_TARGET} src/main_multinode.cpp src/solctra_multinode.cpp src/utils.cpp)
target_include_directories(${EXECUTABLE_TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(MPI REQUIRED)
target_link_libraries(${EXECUTABLE_TARGET} PRIVATE MPI::MPI_CXX)

find_package(Threads REQUIRED)
target_link_libraries(${EXECUTABLE_TARGET} PRIVATE Threads::Threads)

if(${TREAT_WARNINGS_AS_ERRORS})
    target_compile_options(${EXECUTABLE_TARGET} PRIVATE -Werror)
endif()

if(${VECTOR})
    target_compile_options(${EXECUTABLE_TARGET} PRIVATE -ffast-math -march=knl -mavx512f -mavx512pf -mavx512er -mavx512cd)
endif()
