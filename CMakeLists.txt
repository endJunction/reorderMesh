cmake_minimum_required(VERSION 3.1)

project(ReorderMesh)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(Boost REQUIRED COMPONENTS graph)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(reorderMesh main.cpp)

target_link_libraries(reorderMesh ${VTK_LIBRARIES} ${Boost_LIBRARIES})
