cmake_minimum_required (VERSION 3.0)
project (emp-aby)
set(NAME "emp-aby")

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}-O3 -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer -Wno-sign-compare")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenFHE_CXX_FLAGS} -O3 -fno-omit-frame-pointer -mno-omit-leaf-frame-pointer -DMATHBACKEND=4")

if(NOT ${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DLIBS_USE_LIBUV -DLIBS_USE_OPENSSL")
endif()

set(CMAKE_CXX_STANDARD 17)

find_path(CMAKE_FOLDER NAMES cmake/emp-tool-config.cmake)
include(${CMAKE_FOLDER}/cmake/emp-base.cmake)

find_package(OpenFHE REQUIRED)

include_directories(${OPENMP_INCLUDES})
include_directories(${OpenFHE_INCLUDE})
include_directories(${OpenFHE_INCLUDE}/third-party/include)
include_directories(${OpenFHE_INCLUDE}/core)
include_directories(${OpenFHE_INCLUDE}/pke)
link_directories(${OpenFHE_LIBDIR})
link_directories(${OPENMP_LIBRARIES})
if(BUILD_STATIC)
    set( CMAKE_EXE_LINKER_FLAGS "${OpenFHE_EXE_LINKER_FLAGS} -static")
    link_libraries( ${OpenFHE_STATIC_LIBRARIES} )
else()
    set( CMAKE_EXE_LINKER_FLAGS ${OpenFHE_EXE_LINKER_FLAGS} )
    link_libraries( ${OpenFHE_SHARED_LIBRARIES} )
endif()

find_package(emp-ot REQUIRED)
include_directories(${EMP-OT_INCLUDE_DIRS})
# target_link_libraries(emp-aby ${OPENFHE_LIBS})

# Installation
install(FILES cmake/emp-aby-config.cmake DESTINATION cmake/)
install(DIRECTORY emp-aby DESTINATION include/)

ENABLE_TESTING()
ADD_SUBDIRECTORY(test)


