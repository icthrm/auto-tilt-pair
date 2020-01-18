
#set(KEYMATCH_SRCS match_core.cpp)
aux_source_directory(. MRC_IMG_SRCS)

MESSAGE(STATUS "building mrcimg")

add_library (mrcimg STATIC ${MRC_IMG_SRCS})

find_package(OpenCV REQUIRED)

target_link_libraries(mrcimg ${OpenCV_LIBS})#opencv_core opencv_imgproc opencv_calib3d opencv_video opencv_features2d 
#	opencv_ml opencv_highgui opencv_objdetect opencv_contrib opencv_legacy)

set_target_properties(mrcimg PROPERTIES OUTPUT_NAME "mrcimg")
set_target_properties(mrcimg PROPERTIES CLEAN_DIRECT_OUTPUT 1)
