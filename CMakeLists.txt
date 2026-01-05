cmake_minimum_required (VERSION 3.0)
project (emp-aby)
set(NAME "emp-aby")

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}-O3 -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer -Wno-sign-compare")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3 -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer -DMATHBACKEND=4")

set(CMAKE_CXX_STANDARD 17)

find_path(CMAKE_FOLDER NAMES cmake/emp-tool-config.cmake)
include(${CMAKE_FOLDER}/cmake/emp-base.cmake)


include_directories(${OPENMP_INCLUDES})
link_directories(${OPENMP_LIBRARIES})
find_package(emp-ot REQUIRED)
include_directories(${EMP-OT_INCLUDE_DIRS})
# target_link_libraries(emp-aby ${OPENFHE_LIBS})

# Installation
install(FILES cmake/emp-aby-config.cmake DESTINATION cmake/)
install(DIRECTORY emp-aby DESTINATION include/)

ENABLE_TESTING()
ADD_SUBDIRECTORY(test)


