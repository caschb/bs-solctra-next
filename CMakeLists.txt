cmake_minimum_required(VERSION 3.21)

if (NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
endif()

set(CMAKE_CXX_EXTENSIONS OFF)

project(
  BS-Solctra
  VERSION 0.2.0
  DESCRIPTION "biot-savart solver"
  HOMEPAGE_URL "https://gitlab.com/CNCA_CeNAT/bs-solctra-modern"
  LANGUAGES CXX C)

include(cmake/PreventInSourceBuilds.cmake)
include(ProjectOptions.cmake)


bs_solctra_setup_options()

bs_solctra_global_options()
include(Dependencies.cmake)
bs_solctra_setup_dependencies()

bs_solctra_local_options()

# don't know if this should be set globally from here or not...
set(CMAKE_CXX_VISIBILITY_PRESET hidden)

set(GIT_SHA
    "Unknown"
    CACHE STRING "SHA this build was generated from")
string(
  SUBSTRING "${GIT_SHA}"
            0
            8
            GIT_SHORT_SHA)

target_compile_features(bs_solctra_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

add_library(bs_solctra::bs_solctra_options ALIAS bs_solctra_options)
add_library(bs_solctra::bs_solctra_warnings ALIAS bs_solctra_warnings)

# configure files based on CMake configuration options
# add_subdirectory(configured_files)

# Adding the src:
add_subdirectory(main)

# Don't even look at tests if we're not top level
if(NOT PROJECT_IS_TOP_LEVEL)
  return()
endif()

# Adding the tests:
include(CTest)

# if(BUILD_TESTING)
#   message(AUTHOR_WARNING "Building Tests. Be sure to check out test/constexpr_tests.cpp for constexpr testing")
#   add_subdirectory(test)
# endif()


# If MSVC is being used, and ASAN is enabled, we need to set the debugger environment
# so that it behaves well with MSVC's debugger, and we can run the target from visual studio
if(MSVC)
  get_all_installable_targets(all_targets)
  message("all_targets=${all_targets}")
  set_target_properties(${all_targets} PROPERTIES VS_DEBUGGER_ENVIRONMENT "PATH=$(VC_ExecutablePath_x64);%PATH%")
endif()

# set the startup project for the "play" button in MSVC
set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT intro)

if(CMAKE_SKIP_INSTALL_RULES)
  return()
endif()
