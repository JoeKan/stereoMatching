project (StereoMatching)

cmake_minimum_required (VERSION 3.1)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(FREEIMAGE_FIND_REQUIRED, TRUE)
find_package(FreeImage)

add_executable(StereoMatching src/stereo.cpp src/FreeImageHelper.cpp src/Eigen.h src/FreeImageHelper.h)

target_include_directories(StereoMatching PUBLIC "${PROJECT_SOURCE_DIR}/libs/Eigen" ${FREEIMAGE_INCLUDE_DIRS})
target_link_libraries(StereoMatching ${FREEIMAGE_LIBRARIES})