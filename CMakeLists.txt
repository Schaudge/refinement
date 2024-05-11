cmake_minimum_required(VERSION 3.10)
project(recsc)

set(CMAKE_CXX_STANDARD 20)
add_compile_options(-Wall -Wno-unused-function -Wno-unused-variable -Wno-sign-compare)

include_directories(lib)
find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

# discover static libhts, make sure the directory is an absolute path.
set(HTSLIB_DIR ${PROJECT_SOURCE_DIR}/lib)
# set(HTSLIB_DIR /usr/local/lib/) # if htslib is installed system-wide
find_library(HTS_LIB NAMES hts HINTS ${HTSLIB_DIR} STATIC IMPORTED)

# add all source file to the SRC_LIST *variable*
aux_source_directory(src SRC_LIST)

add_executable(recsc src/main.cpp ${SRC_LIST})
target_link_libraries(recsc ${HTS_LIB} z bz2 lzma curl pthread)

