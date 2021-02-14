project(burgercpp)

cmake_minimum_required(VERSION 3.2.0)

# C++14 support
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -std=c++14")

# Sources
set (burgercpp_SRC
        "src/main.cpp"
        "src/mesh.h"
        "src/scene.cpp"
        "src/objloader.h"
        "src/ray.cpp"
        "src/rfimage.h"
        "src/psf.h"
        "src/transducer.h"
        "src/volume.h"
        "src/wavefront/tiny_obj_loader.cpp"
        )

add_executable(burgercpp ${burgercpp_SRC})

# Bullet includes
set(BULLET_LIB_DIR "/usr/local/lib")
set(BULLET_SRC_DIR "/usr/local/include/bullet")
#set(BULLET_SRC_DIR "" CACHE PATH "")
#set(BULLET_LIB_DIR "" CACHE PATH "")
message("Bullet source directory: ${BULLET_SRC_DIR}")
message("Bullet lib directory: ${BULLET_LIB_DIR}")

# EnvelopeExtraction includes
set(FFTW_LIB_DIR "/usr/local/lib")
set(FFTW_SRC_DIR "/usr/local/include")
#set(FFTW_SRC_DIR "" CACHE PATH "")
#set(FFTW_LIB_DIR "" CACHE PATH "")
message("EnvelopeExtraction source directory: ${FFTW_SRC_DIR}")
message("EnvelopeExtraction lib directory: ${FFTW_LIB_DIR}")

include_directories(${BULLET_SRC_DIR})
include_directories(${FFTW_SRC_DIR})

target_link_libraries(burgercpp
    ${BULLET_LIB_DIR}/libBullet3Common.so
    ${BULLET_LIB_DIR}/libBulletDynamics.so
    ${BULLET_LIB_DIR}/libBulletCollision.so
    ${BULLET_LIB_DIR}/libLinearMath.so
    ${FFTW_LIB_DIR}/libfftw3.a
    )

# OpenCV
find_package(OpenCV REQUIRED)
if (OpenCV_FOUND)
  target_link_libraries(burgercpp ${OpenCV_LIBS})
endif()

# Include other stuff, like Units
include_directories("${CMAKE_SOURCE_DIR}/include")
