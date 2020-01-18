
MESSAGE(STATUS "building main")

# find_package(Eigen3 REQUIRED)

# include_directories(${PROJECT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR})

add_executable(drawmatch main.cpp) #cpd_test.cpp)

# list(APPEND CMAKE_CXX_FLAGS "-std=c++0x")

target_link_libraries(drawmatch util dataf mrcimg modelmatch)
