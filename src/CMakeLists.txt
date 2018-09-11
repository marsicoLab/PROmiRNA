cmake_minimum_required (VERSION 3.0.0)

project (PROmiRNA C CXX)

set (SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/)
set (CMAKE_BUILD_TYPE "Release")

# ----------------------------------------------------------------------------
# Requirements and settings for R embedding
# ----------------------------------------------------------------------------

execute_process(COMMAND R RHOME
                OUTPUT_VARIABLE R_HOME)

set(NUM_TRUNC_CHARS 2)

set (RPATH "R")
set (RSCRIPT_PATH "Rscript")

set (RCPPFLAGS_CMD " ${RPATH} " " CMD " " config " " --cppflags ")

execute_process(COMMAND ${RPATH} CMD config --cppflags
                OUTPUT_VARIABLE RCPPFLAGS)

string(SUBSTRING ${RCPPFLAGS} ${NUM_TRUNC_CHARS} -1 RCPPFLAGS)
include_directories(${RCPPFLAGS})

execute_process(COMMAND ${RPATH} CMD config --ldflags
                OUTPUT_VARIABLE RLDFLAGS)
string(LENGTH ${RLDFLAGS} RLDFLAGS_LEN)

if (${RLDFLAGS} MATCHES "[-][L]([^ ;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RLDFLAGS_L)
    string(STRIP ${RLDFLAGS_L} RLDFLAGS_L )
    link_directories(${RLDFLAGS_L} )
endif()

if (${RLDFLAGS} MATCHES "[-][l]([^;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RLDFLAGS_l)
    string(STRIP ${RLDFLAGS_l} RLDFLAGS_l )
endif()

execute_process(COMMAND ${RSCRIPT_PATH} -e "Rcpp:::CxxFlags()"
                OUTPUT_VARIABLE RCPPINCL)
string(SUBSTRING ${RCPPINCL} ${NUM_TRUNC_CHARS} -1 RCPPINCL)
include_directories(${RCPPINCL})

execute_process(COMMAND ${RSCRIPT_PATH} -e "Rcpp:::LdFlags()"
                OUTPUT_VARIABLE RCPPLIBS)

execute_process(COMMAND ${RSCRIPT_PATH} -e "RInside:::CxxFlags()"
                OUTPUT_VARIABLE RINSIDEINCL)
string(SUBSTRING ${RINSIDEINCL} ${NUM_TRUNC_CHARS} -1 RINSIDEINCL)
include_directories(${RINSIDEINCL})
message(${RINSIDEINCL})

execute_process(COMMAND ${RSCRIPT_PATH} -e "RInside:::LdFlags()"
                OUTPUT_VARIABLE RINSIDELIBS)

execute_process(COMMAND ${RSCRIPT_PATH} -e "RcppArmadillo:::CxxFlags()"
                OUTPUT_VARIABLE RCPPARMADILLOINCL)

string(LENGTH ${RCPPARMADILLOINCL} INCLLENGTH)
math(EXPR INCLLENGTH "${INCLLENGTH}-4")

string(SUBSTRING ${RCPPARMADILLOINCL} 3 ${INCLLENGTH} RCPPARMADILLOINCL)
include_directories(${RCPPARMADILLOINCL})

if (${RCPPLIBS} MATCHES "[-][L]([^ ;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RCPPLIBS_L)
    link_directories(${RCPPLIBS_L} )
endif()

if (${RCPPLIBS} MATCHES "[-][l][R]([^;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RCPPLIBS_l)
endif()

if (${RINSIDELIBS} MATCHES "[-][L]([^ ;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RINSIDELIBS_L)
    link_directories(${RINSIDELIBS_L})
endif()

if (${RINSIDELIBS} MATCHES "[-][l][R]([^;])+")
    string(SUBSTRING ${CMAKE_MATCH_0} ${NUM_TRUNC_CHARS} -1 RINSIDELIBS_l)
endif()

# If the faster 'gold' linker is used, to avoid complaints about undefined symbol
SET(CMAKE_FIND_LIBRARY_SUFFIXES_SAVED ${CMAKE_FIND_LIBRARY_SUFFIXES}) # Backup
LIST(APPEND CMAKE_FIND_LIBRARY_SUFFIXES ".so.3")
FIND_LIBRARY(BLAS_LIBRARY blas)
SET(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAVED}) # Restore
LIST(APPEND EXTRA_LIBRARIES ${BLAS_LIBRARY})

FIND_LIBRARY(LAPACK_LIBRARY lapack)

execute_process(COMMAND ${RPATH} CMD config CXXFLAGS
                OUTPUT_VARIABLE RCXXFLAGS)

execute_process(COMMAND ${RPATH} CMD config BLAS_LIBS
                OUTPUT_VARIABLE RBLAS)

execute_process(COMMAND ${RPATH} CMD config LAPACK_LIBS
                OUTPUT_VARIABLE RLAPACK)

# ----------------------------------------------------------------------------
# SeqAn and TRAP settings
# ----------------------------------------------------------------------------

set (CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../seqan/util/cmake")
set (SEQAN_INCLUDE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../seqan/include")
set (TRAP_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/TRAP/)

find_package (OpenMP REQUIRED)
find_package (SeqAn REQUIRED)
find_package (ZLIB REQUIRED)

include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})

# ----------------------------------------------------------------------------
# Old includes and make calls for hgWiggle from kentUtils
# ----------------------------------------------------------------------------

#set (kentUtils_path ${CMAKE_CURRENT_SOURCE_DIR}/../kentUtils/)

#add_custom_target(build_kentutils_fetch_samtabix make "installEnvironment" "fetchSamtabix" WORKING_DIRECTORY ${kentUtils_path})

#add_custom_target(build_kentutils_jkweb_lib make WORKING_DIRECTORY ${kentUtils_path}/src/lib)
#add_dependencies(build_kentutils_jkweb_lib build_kentutils_fetch_samtabix)

#add_custom_target(build_kentutils_jkhgap_lib make WORKING_DIRECTORY ${kentUtils_path}/src/hg/lib)
#add_dependencies(build_kentutils_jkhgap_lib build_kentutils_fetch_samtabix)
#add_dependencies(build_kentutils_jkhgap_lib build_kentutils_jkweb_lib)

#add_custom_target(build_kentutils_samtabix_lib make libsamtabix.a WORKING_DIRECTORY ${kentUtils_path}/samtabix)
#add_dependencies(build_kentutils_samtabix_lib build_kentutils_fetch_samtabix)

#add_library(kentutils_jkweb_lib STATIC IMPORTED)
#set_target_properties(kentutils_jkweb_lib PROPERTIES IMPORTED_LOCATION ${kentUtils_path}/src/lib/local/jkweb.a)
#add_dependencies(kentutils_jkweb_lib build_kentutils_jkweb_lib)

#add_library(kentutils_jkhgap_lib STATIC IMPORTED)
#set_target_properties(kentutils_jkhgap_lib PROPERTIES IMPORTED_LOCATION ${kentUtils_path}/src/lib/local/jkhgap.a)
#add_dependencies(kentutils_jkhgap_lib build_kentutils_jkhgap_lib)

#add_library(kentutils_samtabix_lib STATIC IMPORTED)
#set_target_properties(kentutils_samtabix_lib PROPERTIES IMPORTED_LOCATION ${kentUtils_path}/samtabix/libsamtabix.a)
#add_dependencies(kentutils_samtabix_lib build_kentutils_samtabix_lib)

#add_library(hgWiggle ${CMAKE_CURRENT_SOURCE_DIR}/hgWiggle/hgWiggleAdapted.c)
#target_link_libraries(hgWiggle kentutils_jkhgap_lib kentutils_jkweb_lib pthread ssl crypto kentutils_samtabix_lib z m mysqlclient)
#target_include_directories(hgWiggle PUBLIC ${kentUtils_path}/src/inc ${kentUtils_path}/src/hg/inc)

# ----------------------------------------------------------------------------
# PROmiRNA executable including CD-HIT and TRAP source files
# ----------------------------------------------------------------------------

set (CMAKE_CXX_FLAGS "-std=c++17 -O3 ${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

message ("cmake_cxx_flags" ${CMAKE_CXX_FLAGS})

add_executable (PROmiRNA PROmiRNA.cpp ${CMAKE_CURRENT_SOURCE_DIR}/CD-HIT/cdhit-common.c++ ${CMAKE_CURRENT_SOURCE_DIR}/CD-HIT/cdhit-utility.c++ ${CMAKE_CURRENT_SOURCE_DIR}/CD-HIT/cdhit-est.c++ ${CMAKE_CURRENT_SOURCE_DIR}/TRAP/annotate_sequence.cpp ${CMAKE_CURRENT_SOURCE_DIR}/TRAP/matrix_properties.cpp ${CMAKE_CURRENT_SOURCE_DIR}/TRAP/read_matrix.cpp ${CMAKE_CURRENT_SOURCE_DIR}/TRAP/region_read.cpp ${CMAKE_CURRENT_SOURCE_DIR}/TRAP/score_calc.cpp)
target_link_libraries (PROmiRNA ${SEQAN_LIBRARIES} "stdc++fs")
#target_link_libraries (PROmiRNA hgWiggle)
target_link_libraries (PROmiRNA)
# Libraries needed for R embedding
target_link_libraries(PROmiRNA ${RLDFLAGS_l})
target_link_libraries(PROmiRNA ${BLAS_LIBS})
target_link_libraries(PROmiRNA ${LAPACK_LIBS})
target_link_libraries(PROmiRNA ${RINSIDELIBS_l})
target_link_libraries(PROmiRNA ${RCPPLIBS_l})
target_link_libraries(PROmiRNA ${EXTRA_LIBRARIES})
target_link_libraries(PROmiRNA ${LAPACK_LIBRARY})
# Include directory for TRAP header files
target_include_directories (PROmiRNA PUBLIC ${TRAP_INCLUDE_PATH})
#add_dependencies(PROmiRNA hgWiggle)

if(OLD_MIRBASE_VERSION)
    target_compile_definitions(PROmiRNA PUBLIC OLD_MIRBASE_VERSION=1)
endif()
