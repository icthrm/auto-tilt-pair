#set(DATAF_SRCS )
aux_source_directory(. DATAF_SRCS)

MESSAGE(STATUS "building dataf")

include_directories(${PROJECT_SOURCE_DIR})

add_library (dataf STATIC ${DATAF_SRCS})

target_link_libraries(dataf util)

set_target_properties(dataf PROPERTIES OUTPUT_NAME "dataf")
set_target_properties(dataf PROPERTIES CLEAN_DIRECT_OUTPUT 1)
