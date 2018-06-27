project (StereoMatching)

cmake_minimum_required (VERSION 3.1)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(FREEIMAGE_FIND_REQUIRED, TRUE)
set(OPENEXR_INCLUDE, "/media/virendra/data/study/4_sem/3D_Scan/Project/stereoMatching/libs/openexr/include")
set(OPENEXR_LIB, "/media/virendra/data/study/4_sem/3D_Scan/Project/stereoMatching/libs/openexr/lib")
find_package(FreeImage)
find_package(TBB)

add_executable(StereoMatching src/stereo.cpp src/FreeImageHelper.cpp src/Eigen.h src/FreeImageHelper.h)

target_include_directories(StereoMatching PRIVATE "${PROJECT_SOURCE_DIR}/libs/Eigen" ${FREEIMAGE_INCLUDE_DIRS} ${TBB_INCLUDE_DIRS} ${OPENEXR_INCLUDE})
target_link_libraries(StereoMatching ${FREEIMAGE_LIBRARIES} ${TBB_LIBRARIES} ${OPENEXR_LIB})