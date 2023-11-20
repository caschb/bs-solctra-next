set(EXECUTABLE_TARGET bs_solctra)

if(bs_solctra_USE_ASCENT)
    add_executable(${EXECUTABLE_TARGET} src/main_multinode.cpp src/solctra_multinode.cpp src/utils.cpp src/ascent_adaptor.cpp)
elseif(bs_solctra_USE_CATALYST)
    add_executable(${EXECUTABLE_TARGET} src/main_multinode.cpp src/solctra_multinode.cpp src/utils.cpp src/catalyst_adaptor.cpp)
else()
    add_executable(${EXECUTABLE_TARGET} src/main_multinode.cpp src/solctra_multinode.cpp src/utils.cpp)
endif()


target_link_libraries(
  ${EXECUTABLE_TARGET}
  PRIVATE bs_solctra::bs_solctra_options
          bs_solctra::bs_solctra_warnings)


target_include_directories(${EXECUTABLE_TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(MPI REQUIRED)
target_link_system_libraries(${EXECUTABLE_TARGET} PRIVATE MPI::MPI_CXX)

find_package(OpenMP REQUIRED)
target_link_system_libraries(${EXECUTABLE_TARGET} PRIVATE OpenMP::OpenMP_CXX)

find_package(spdlog REQUIRED)
target_link_system_libraries(${EXECUTABLE_TARGET} PRIVATE spdlog::spdlog)

find_package(CLI11 REQUIRED)
target_link_system_libraries(${EXECUTABLE_TARGET} PRIVATE CLI11::CLI11)

if(bs_solctra_USE_ASCENT)
    find_package(Ascent REQUIRED
        NO_DEFAULT_PATH
        PATHS "/home/casch/Devel/ascent/install-debug")
    target_link_system_libraries(${EXECUTABLE_TARGET} PRIVATE ascent::ascent)
elseif(bs_solctra_USE_CATALYST)
    find_package(catalyst 2.0 REQUIRED)
    target_link_system_libraries(${EXECUTABLE_TARGET} PRIVATE catalyst::catalyst)
endif()

if(bs_solctra_solctra_ENABLE_KNL_OPTS)
    target_compile_options(${EXECUTABLE_TARGET} PRIVATE -mavx512f -mavx512pf -mavx512er -mavx512cd -march=knl)
endif()
if(bs_solctra_ENABLE_FAST_MATH)
    target_compile_options(${EXECUTABLE_TARGET} PRIVATE -ffast-math)
endif()

# if(${GENERATE_OPTIMIZATION_REPORTS})
#     if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
#     target_compile_options(${EXECUTABLE_TARGET} PRIVATE -Rpass="*" -Rpass-missed="*" -Rpass-analysis="*" -fsave-optimization-record=yaml -foptimization-record-file=report-clang)
#     endif()

#     if(${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
#     target_compile_options(${EXECUTABLE_TARGET} PRIVATE -fopt-info-optall=report)
#     endif()
# endif()