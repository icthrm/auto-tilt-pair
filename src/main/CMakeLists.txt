
MESSAGE(STATUS "building main")

# include_directories(${PROJECT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR})

add_executable(autopairs main.cpp) #cpd_test.cpp)

# list(APPEND CMAKE_CXX_FLAGS "-std=c++0x")

target_link_libraries(autopairs util dataf modelmatch)
