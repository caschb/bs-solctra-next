set(EXECUTABLE_TARGET bs_solctra)

add_executable(${EXECUTABLE_TARGET} src/main_multinode.cpp src/solctra_multinode.cpp src/utils.cpp)
target_include_directories(${EXECUTABLE_TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(MPI REQUIRED)
target_link_libraries(${EXECUTABLE_TARGET} PRIVATE MPI::MPI_CXX)

find_package(OpenMP REQUIRED)
target_link_libraries(${EXECUTABLE_TARGET} PRIVATE OpenMP::OpenMP_CXX)

if(${TREAT_WARNINGS_AS_ERRORS})
    target_compile_options(${EXECUTABLE_TARGET} PRIVATE -Werror)
endif()

if(${VECTOR})
    target_compile_options(${EXECUTABLE_TARGET} PRIVATE -mavx512f -mavx512pf -mavx512er -mavx512cd)
endif()

if(${FAST_MATH})
    target_compile_options(${EXECUTABLE_TARGET} PRIVATE -ffast-math)
endif()

if(${ACTIVATE_SANITIZERS})
    target_compile_options(${EXECUTABLE_TARGET} PRIVATE -fsanitize=address)
    target_link_options(${EXECUTABLE_TARGET} PRIVATE -fsanitize=address)
endif()

if(${GENERATE_OPTIMIZATION_REPORTS})
    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    target_compile_options(${EXECUTABLE_TARGET} PRIVATE -Rpass="*" -Rpass-missed="*" -Rpass-analysis="*" -fsave-optimization-record=yaml -foptimization-record-file=report-clang)
    endif()

    if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(${EXECUTABLE_TARGET} PRIVATE -fopt-info-optall=report)
    endif()
endif()