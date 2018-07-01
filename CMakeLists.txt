project (StereoMatching)

cmake_minimum_required (VERSION 3.1)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(FREEIMAGE_FIND_REQUIRED, TRUE)
find_package(FreeImage)
find_package(TBB)
#find_package(llmBase)
find_package(OpenEXR)

add_executable(StereoMatching src/stereo.cpp src/FreeImageHelper.cpp src/Eigen.h src/FreeImageHelper.h src/BruteForce.cpp src/parameters.h src/utility.h src/utility.cpp)

target_include_directories(StereoMatching PRIVATE "${PROJECT_SOURCE_DIR}/libs/Eigen" ${FREEIMAGE_INCLUDE_DIRS} ${TBB_INCLUDE_DIRS} ${OPENEXR_INCLUDE_DIR})
target_link_libraries(StereoMatching ${FREEIMAGE_LIBRARIES} ${TBB_LIBRARIES} ${OPENEXR_LIBRARIES})