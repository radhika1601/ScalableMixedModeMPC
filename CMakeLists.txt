cmake_minimum_required(VERSION 3.25)
project(emp-aby LANGUAGES C CXX)
set(NAME "emp-aby")

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -Wno-sign-compare")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Apple clang searches /usr/local/include before -isystem paths, so keep emp includes as -I.
set(CMAKE_NO_SYSTEM_FROM_IMPORTED ON)
find_package(emp-ot CONFIG REQUIRED)

include_directories(${CMAKE_SOURCE_DIR})

add_library(emp-aby INTERFACE)
target_include_directories(emp-aby INTERFACE
	$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include>)
target_link_libraries(emp-aby INTERFACE emp-ot::emp-ot)

# Installation
install(FILES cmake/emp-aby-config.cmake DESTINATION cmake/)
install(DIRECTORY emp-aby DESTINATION include/)

ENABLE_TESTING()
ADD_SUBDIRECTORY(test)
