include(cmake/SystemLink.cmake)
include(CMakeDependentOption)
include(CheckCXXCompilerFlag)


macro(bs_solctra_supports_sanitizers)
  if((CMAKE_CXX_COMPILER_ID MATCHES ".*Clang.*" OR CMAKE_CXX_COMPILER_ID MATCHES ".*GNU.*") AND NOT WIN32)
    set(SUPPORTS_UBSAN ON)
  else()
    set(SUPPORTS_UBSAN OFF)
  endif()

  if((CMAKE_CXX_COMPILER_ID MATCHES ".*Clang.*" OR CMAKE_CXX_COMPILER_ID MATCHES ".*GNU.*") AND WIN32)
    set(SUPPORTS_ASAN OFF)
  else()
    set(SUPPORTS_ASAN ON)
  endif()
endmacro()

macro(bs_solctra_setup_options)
  option(bs_solctra_ENABLE_HARDENING "Enable hardening" ON)
  option(bs_solctra_ENABLE_COVERAGE "Enable coverage reporting" OFF)
  cmake_dependent_option(
    bs_solctra_ENABLE_GLOBAL_HARDENING
    "Attempt to push hardening options to built dependencies"
    ON
    bs_solctra_ENABLE_HARDENING
    OFF)

  bs_solctra_supports_sanitizers()

  if(NOT PROJECT_IS_TOP_LEVEL OR bs_solctra_PACKAGING_MAINTAINER_MODE)
    option(bs_solctra_ENABLE_IPO "Enable IPO/LTO" OFF)
    option(bs_solctra_WARNINGS_AS_ERRORS "Treat Warnings As Errors" OFF)
    option(bs_solctra_ENABLE_USER_LINKER "Enable user-selected linker" OFF)
    option(bs_solctra_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" OFF)
    option(bs_solctra_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
    option(bs_solctra_ENABLE_SANITIZER_UNDEFINED "Enable undefined sanitizer" OFF)
    option(bs_solctra_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
    option(bs_solctra_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)
    option(bs_solctra_ENABLE_UNITY_BUILD "Enable unity builds" OFF)
    option(bs_solctra_ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)
    option(bs_solctra_ENABLE_CPPCHECK "Enable cpp-check analysis" OFF)
    option(bs_solctra_ENABLE_PCH "Enable precompiled headers" OFF)
    option(bs_solctra_ENABLE_CACHE "Enable ccache" OFF)
    option(bs_solctra_ENABLE_KNL_OPTS "Enable KNL optimizations" OFF)
    option(bs_solctra_ENABLE_FAST_MATH "Enable fast math" OFF)
    option(bs_solctra_USE_CATALYST "Use Catalyst" OFF)
    option(bs_solctra_USE_ASCENT "Use Ascent" OFF)
  else()
    option(bs_solctra_ENABLE_IPO "Enable IPO/LTO" ON)
    option(bs_solctra_WARNINGS_AS_ERRORS "Treat Warnings As Errors" ON)
    option(bs_solctra_ENABLE_USER_LINKER "Enable user-selected linker" OFF)
    option(bs_solctra_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" ${SUPPORTS_ASAN})
    option(bs_solctra_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
    option(bs_solctra_ENABLE_SANITIZER_UNDEFINED "Enable undefined sanitizer" ${SUPPORTS_UBSAN})
    option(bs_solctra_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
    option(bs_solctra_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)
    option(bs_solctra_ENABLE_UNITY_BUILD "Enable unity builds" OFF)
    option(bs_solctra_ENABLE_CLANG_TIDY "Enable clang-tidy" ON)
    option(bs_solctra_ENABLE_CPPCHECK "Enable cpp-check analysis" ON)
    option(bs_solctra_ENABLE_PCH "Enable precompiled headers" OFF)
    option(bs_solctra_ENABLE_CACHE "Enable ccache" ON)
    option(bs_solctra_ENABLE_KNL_OPTS "Enable KNL optimizations" OFF)
    option(bs_solctra_ENABLE_FAST_MATH "Enable fast math" OFF)
    option(bs_solctra_USE_CATALYST "Use Catalyst" OFF)
    option(bs_solctra_USE_ASCENT "Use Ascent" OFF)
  endif()

  if(NOT PROJECT_IS_TOP_LEVEL)
    mark_as_advanced(
      bs_solctra_ENABLE_IPO
      bs_solctra_WARNINGS_AS_ERRORS
      bs_solctra_ENABLE_USER_LINKER
      bs_solctra_ENABLE_SANITIZER_ADDRESS
      bs_solctra_ENABLE_SANITIZER_LEAK
      bs_solctra_ENABLE_SANITIZER_UNDEFINED
      bs_solctra_ENABLE_SANITIZER_THREAD
      bs_solctra_ENABLE_SANITIZER_MEMORY
      bs_solctra_ENABLE_UNITY_BUILD
      bs_solctra_ENABLE_CLANG_TIDY
      bs_solctra_ENABLE_CPPCHECK
      bs_solctra_ENABLE_COVERAGE
      bs_solctra_ENABLE_PCH
      bs_solctra_ENABLE_CACHE)
  endif()
endmacro()

macro(bs_solctra_global_options)
  if(bs_solctra_ENABLE_IPO)
    include(cmake/InterproceduralOptimization.cmake)
    bs_solctra_enable_ipo()
  endif()

  bs_solctra_supports_sanitizers()

  if(bs_solctra_ENABLE_HARDENING AND bs_solctra_ENABLE_GLOBAL_HARDENING)
    include(cmake/Hardening.cmake)
    if(NOT SUPPORTS_UBSAN 
       OR bs_solctra_ENABLE_SANITIZER_UNDEFINED
       OR bs_solctra_ENABLE_SANITIZER_ADDRESS
       OR bs_solctra_ENABLE_SANITIZER_THREAD
       OR bs_solctra_ENABLE_SANITIZER_LEAK)
      set(ENABLE_UBSAN_MINIMAL_RUNTIME FALSE)
    else()
      set(ENABLE_UBSAN_MINIMAL_RUNTIME TRUE)
    endif()
    bs_solctra_enable_hardening(bs_solctra_options ON ${ENABLE_UBSAN_MINIMAL_RUNTIME})
  endif()
endmacro()

macro(bs_solctra_local_options)
  if(PROJECT_IS_TOP_LEVEL)
    include(cmake/StandardProjectSettings.cmake)
  endif()

  add_library(bs_solctra_warnings INTERFACE)
  add_library(bs_solctra_options INTERFACE)

  include(cmake/CompilerWarnings.cmake)
  bs_solctra_set_project_warnings(
    bs_solctra_warnings
    ${bs_solctra_WARNINGS_AS_ERRORS}
    ""
    ""
    ""
    "")

  if(bs_solctra_ENABLE_USER_LINKER)
    include(cmake/Linker.cmake)
    configure_linker(bs_solctra_options)
  endif()

  include(cmake/Sanitizers.cmake)
  bs_solctra_enable_sanitizers(
    bs_solctra_options
    ${bs_solctra_ENABLE_SANITIZER_ADDRESS}
    ${bs_solctra_ENABLE_SANITIZER_LEAK}
    ${bs_solctra_ENABLE_SANITIZER_UNDEFINED}
    ${bs_solctra_ENABLE_SANITIZER_THREAD}
    ${bs_solctra_ENABLE_SANITIZER_MEMORY})

  set_target_properties(bs_solctra_options PROPERTIES UNITY_BUILD ${bs_solctra_ENABLE_UNITY_BUILD})

  if(bs_solctra_ENABLE_PCH)
    target_precompile_headers(
      bs_solctra_options
      INTERFACE
      <vector>
      <string>
      <utility>)
  endif()

  if(bs_solctra_ENABLE_CACHE)
    include(cmake/Cache.cmake)
    bs_solctra_enable_cache()
  endif()

  include(cmake/StaticAnalyzers.cmake)
  if(bs_solctra_ENABLE_CLANG_TIDY)
    bs_solctra_enable_clang_tidy(bs_solctra_options ${bs_solctra_WARNINGS_AS_ERRORS})
  endif()

  if(bs_solctra_ENABLE_CPPCHECK)
    bs_solctra_enable_cppcheck(${bs_solctra_WARNINGS_AS_ERRORS} "" # override cppcheck options
    )
  endif()

  if(bs_solctra_ENABLE_COVERAGE)
    include(cmake/Tests.cmake)
    bs_solctra_enable_coverage(bs_solctra_options)
  endif()

  if(bs_solctra_WARNINGS_AS_ERRORS)
    check_cxx_compiler_flag("-Wl,--fatal-warnings" LINKER_FATAL_WARNINGS)
    if(LINKER_FATAL_WARNINGS)
      # This is not working consistently, so disabling for now
      # target_link_options(bs_solctra_options INTERFACE -Wl,--fatal-warnings)
    endif()
  endif()

  if(bs_solctra_ENABLE_HARDENING AND NOT bs_solctra_ENABLE_GLOBAL_HARDENING)
    include(cmake/Hardening.cmake)
    if(NOT SUPPORTS_UBSAN
       OR bs_solctra_ENABLE_SANITIZER_UNDEFINED
       OR bs_solctra_ENABLE_SANITIZER_ADDRESS
       OR bs_solctra_ENABLE_SANITIZER_THREAD
       OR bs_solctra_ENABLE_SANITIZER_LEAK)
      set(ENABLE_UBSAN_MINIMAL_RUNTIME FALSE)
    else()
      set(ENABLE_UBSAN_MINIMAL_RUNTIME TRUE)
    endif()
    bs_solctra_enable_hardening(bs_solctra_options OFF ${ENABLE_UBSAN_MINIMAL_RUNTIME})
  endif()
endmacro()