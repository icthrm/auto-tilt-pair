aux_source_directory(. MODELMATCH_SRCS)

MESSAGE(STATUS "building modelmatch")

add_library (modelmatch STATIC ${MODELMATCH_SRCS})

link_directories (${PROJECT_BINARY_DIR}/lib)

target_link_libraries(modelmatch blaswrap f2c clapack m util ann_1.1_char mrcimg matrix levmar)

set_target_properties(modelmatch PROPERTIES OUTPUT_NAME "modelmatch")
set_target_properties(modelmatch PROPERTIES CLEAN_DIRECT_OUTPUT 1)