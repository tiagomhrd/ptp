cmake_minimum_required(VERSION 3.23.1)

project(ptp)

add_subdirectory(ptp)
add_subdirectory(third_party/Catch2)
add_subdirectory(tests)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT tests)