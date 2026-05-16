
cmake_minimum_required(VERSION 3.14)

project(THERMINATOR CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(ROOT REQUIRED CONFIG
  HINTS $ENV{ROOTSYS}/cmake
        /opt/local/libexec/root6/share/root/cmake)

include(${ROOT_USE_FILE})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-reorder")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_subdirectory(src)
add_subdirectory(work)

