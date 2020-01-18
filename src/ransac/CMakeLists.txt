
aux_source_directory(. RANSAC_SRCS)

MESSAGE(STATUS "building ransac")

add_library (ransac STATIC ${RANSAC_SRCS})

find_package(OpenCV REQUIRED)

target_link_libraries(ransac ${OpenCV_LIBS})#opencv_core opencv_imgproc opencv_calib3d opencv_video opencv_features2d 
#	opencv_ml opencv_highgui opencv_objdetect opencv_contrib opencv_legacy)

set_target_properties(ransac PROPERTIES OUTPUT_NAME "ransac")
set_target_properties(ransac PROPERTIES CLEAN_DIRECT_OUTPUT 1)

