cmake_minimum_required (VERSION 3.2)
project (PROmiRNA_tests CXX)

# ----------------------------------------------------------------------------
# Use CCache if possible
# ----------------------------------------------------------------------------

find_program (CCACHE_PROGRAM ccache)
if (CCACHE_PROGRAM)
    set_property (GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CCACHE_PROGRAM}")
    set_property (GLOBAL PROPERTY RULE_LAUNCH_LINK "${CCACHE_PROGRAM}")
endif (CCACHE_PROGRAM)

# ----------------------------------------------------------------------------
# Add GTest
# ----------------------------------------------------------------------------

enable_testing ()

set (GOOGLETEST_ROOT "${PROJECT_BINARY_DIR}/vendor/googletest")
include_directories ("${GOOGLETEST_ROOT}/googletest/include/")

include (ExternalProject)
ExternalProject_Add(
    googletest
    PREFIX googletest
    GIT_REPOSITORY "https://github.com/google/googletest.git"
    GIT_TAG "15392f1a38fa0b8c3f13a9732e94b209069efa1c"
    SOURCE_DIR "${GOOGLETEST_ROOT}"
    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR};-DCMAKE_INSTALL_LIBDIR=${PROJECT_BINARY_DIR}/lib/;-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER};-DBUILD_GTEST=1;-DBUILD_GMOCK=0"
    UPDATE_DISCONNECTED yes
)

add_library (gtest_main STATIC IMPORTED)
add_dependencies (gtest_main googletest)
set_target_properties (gtest_main PROPERTIES IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/lib/libgtest_main.a")

add_library (gtest STATIC IMPORTED)
add_dependencies (gtest gtest_main)
set_target_properties (gtest PROPERTIES IMPORTED_LOCATION "${PROJECT_BINARY_DIR}/lib/libgtest.a")

# ----------------------------------------------------------------------------
# Add SeqAn and TRAP
# ----------------------------------------------------------------------------

set (CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../seqan/util/cmake")
set (SEQAN_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../seqan/include")

find_package (OpenMP REQUIRED)
find_package (SeqAn REQUIRED)
find_package (ZLIB REQUIRED)

include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})

set (CMAKE_CXX_FLAGS "-std=c++17 -O3 ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# ----------------------------------------------------------------------------
# Add the tests
# ----------------------------------------------------------------------------

file(GLOB TEST_SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

foreach(_test_file ${TEST_SRC_FILES})
    get_filename_component(_test_name ${_test_file} NAME_WE)
    add_executable(${_test_name} ${_test_file})
    target_link_libraries(${_test_name} gtest gtest_main ${SEQAN_LIBRARIES} "stdc++fs")
    add_test(${_test_name} ${_test_name})
    if(TESTCASE)
        target_compile_definitions(${_test_name} PUBLIC TESTCASE=1)
    endif()
endforeach()
