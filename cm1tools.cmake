LIST(APPEND CM1TOOLS_SRC_FILES
    cm1tools-3.0/parsedir.c
    cm1tools-3.0/hdfio.c
    cm1tools-3.0/readmult.c
    )
SET_SOURCE_FILES_PROPERTIES( ${CM1TOOLS_SRC_FILES} PROPERTIES 
    COMPILE_DEFINITIONS UNIX
    LANGUAGE C
    COMPILE_FLAGS "-I${CMAKE_CURRENT_SOURCE_DIR}/cm1tools-3.0")


#
# Augment the LIBM_SOURCES and LIBE_SOURCES with the additional files that
# we need to include in those libraries.
#    
SET(LIBM_SOURCES 
${LIBM_SOURCES}
${CM1TOOLS_SRC_FILES}
)

SET(LIBE_SOURCES 
${LIBE_SOURCES}
${CM1TOOLS_SRC_FILES}
)
