cmake_minimum_required(VERSION 4.0)

project(ReorderMesh CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    IOXML
)

find_package(Boost REQUIRED COMPONENTS graph)

add_executable(reorderMesh main.cpp)

target_compile_features(reorderMesh PRIVATE cxx_std_23)

target_link_libraries(reorderMesh PRIVATE ${VTK_LIBRARIES} Boost::graph)

vtk_module_autoinit(
    TARGETS reorderMesh
    MODULES ${VTK_LIBRARIES}
)
