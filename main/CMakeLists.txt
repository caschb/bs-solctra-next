find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)
add_executable(bs_solctra src/main_multinode.cpp src/solctra_multinode.cpp src/utils.cpp)
target_include_directories(bs_solctra PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(bs_solctra PRIVATE ${OpenMP_CXX_DIRECTORIES})
target_include_directories(bs_solctra PRIVATE ${MPI_CXX_INCLUDE_DIRECTORIES})
target_link_libraries(bs_solctra PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(bs_solctra PRIVATE MPI::MPI_CXX)

set_target_properties(bs_solctra PROPERTIES 
    CXX_EXTENSIONS OFF 
    CXX_STANDARD 11)

if(${TREAT_WARNINGS_AS_ERRORS})
    target_compile_options(bs_solctra PRIVATE -Werror)
endif()

if(${VECTOR})
    target_compile_options(bs_solctra PRIVATE -ffast-math -march=knl -mavx512f -mavx512pf -mavx512er -mavx512cd)
endif()