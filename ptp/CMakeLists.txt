cmake_minimum_required(VERSION 3.23.1)

project(ptp)

add_library(${PROJECT_NAME} STATIC
${CMAKE_CURRENT_SOURCE_DIR}/src/ptp.h 
${CMAKE_CURRENT_SOURCE_DIR}/src/ptp.cpp 
)

set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)

target_include_directories(${PROJECT_NAME} 
PRIVATE src ../third_party/ 
)