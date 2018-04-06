project(DBoW2 CXX)

cmake_minimum_required(VERSION 2.7)
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries?")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -mtune=native -march=native")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

find_package(OpenCV 3 REQUIRED)

set(LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib")

include_directories(include include/DBoW2 include/DUtils ${OpenCV_INCLUDE_DIRS})

set(SRC_FILES
  src/DBoW2/BowVector.cpp
  src/DBoW2/FeatureVector.cpp
  src/DBoW2/FORB.cpp
  src/DBoW2/QueryResults.cpp
  src/DBoW2/ScoringObject.cpp

  src/DUtils/Random.cpp
  src/DUtils/Timestamp.cpp
)

add_library(DBoW2 SHARED ${SRC_FILES})
target_link_libraries(DBoW2 ${OpenCV_LIBS})
