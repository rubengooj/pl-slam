project(PL-SLAM)
cmake_minimum_required(VERSION 2.7)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules")

add_subdirectory(3rdparty/line_descriptor)
add_subdirectory(3rdparty/DBoW2)

set(HAS_MRPT ON CACHE BOOL "Build the PointGrey Bumblebee2 SVO application that employs the MRPT library")

find_package(OpenCV 3 REQUIRED)
find_package(Boost    REQUIRED COMPONENTS regex thread system filesystem)
find_package(g2o      REQUIRED)
find_package(Eigen3   REQUIRED)
find_package(yaml-cpp REQUIRED CONFIG PATHS ${YAML_PATHS}) # YAML library

# MRPT library (optional, only with representation purposes)
if(HAS_MRPT)
    find_package(MRPT 1.9.9 REQUIRED opengl gui hwdrivers)
    add_definitions(-DHAS_MRPT)
endif()

if(COMMAND cmake_policy)
  	cmake_policy(SET CMP0003 NEW)	
endif()
link_directories(${OpenCV_LIBS_DIR})
include_directories(${OpenCV2_INCLUDE_DIRS})

# Odometry library (Future TODO: integrate with this library)
SET( StVO_LIBRARY 		"" CACHE FILEPATH "Visual Odometry Library" )
SET( StVO_INCLUDE_DIR 	"" CACHE PATH "Visual Odometry Include folder" )	


SET(BUILD_SHARED_LIBS ON)
SET(CMAKE_MODULE_PATH $ENV{CMAKE_MODULE_PATH})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -mtune=native -march=native")

add_definitions(-DBOOST_NO_CXX11_SCOPED_ENUMS)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

SET(G2O_LIBS
    g2o::core
    g2o::types_slam2d
    g2o::types_slam3d
    g2o::solver_cholmod
    )

# Include dirs
include_directories(
  include
  ${Eigen3_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${YAML_CPP_INCLUDE_DIRS}
  ${PROJECT_SOURCE_DIR}/../stvo-pl/include/
  ${PROJECT_SOURCE_DIR}/3rdparty/DBoW2/include/
  ${PROJECT_SOURCE_DIR}/3rdparty/line_descriptor/include/
  ${G2O_INCLUDE_DIR}
  ${g2o_INCLUDE_DIR}
  /usr/include/suitesparse # for cholmod
)


# Set link libraries
list(APPEND LINK_LIBS 
  ${OpenCV_LIBS}
  ${Boost_LIBRARIES}
  ${YAML_CPP_LIBRARIES}
  ${G2O_LIBS}
  ${PROJECT_SOURCE_DIR}/../stvo-pl/lib/libstvo.so
  DBoW2
  linedesc
)

list(APPEND SOURCEFILES
  src/mapHandler.cpp
  src/mapFeatures.cpp
  src/keyFrame.cpp
  src/slamConfig.cpp
)

# Set source files 
if(HAS_MRPT)
    list(APPEND SOURCEFILES
        src/slamScene.cpp
    )
endif()

# List all files (headers) contained by StVO-PL library
file(GLOB_RECURSE all_include_files RELATIVE "${CMAKE_SOURCE_DIR}" *.h *.hpp)

# Visualize the files of this directory in IDE creating an custom empty target
add_custom_target( plslam_includes DEPENDS ${all_include_files} SOURCES ${all_include_files} )

# Create StVO-PL library
add_library(plslam SHARED ${SOURCEFILES})

target_link_libraries(plslam ${LINK_LIBS})

if(HAS_MRPT)
    target_link_libraries(plslam ${MRPT_LIBS} ${MRPT_LIBRARIES})
endif()

# Applications 
if(HAS_MRPT)
    add_executable       ( plslam_dataset app/plslam_dataset.cpp )
    target_link_libraries( plslam_dataset plslam )
endif()







