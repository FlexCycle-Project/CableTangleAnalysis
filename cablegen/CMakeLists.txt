set(KYMOKNOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../KymoKnot-clone)
set(KYMOKNOT_QHULL_DIR ${KYMOKNOT_DIR}/qhull-2012.1)

add_library(kymoknot_qhull STATIC
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/global.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/stat.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/geom2.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/poly2.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/merge.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/libqhull.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/geom.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/poly.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/qset.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/mem.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/random.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/usermem.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/userprintf.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/io.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/user.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/rboxlib.c
    ${KYMOKNOT_QHULL_DIR}/src/libqhull/userprintf_rbox.c
)
target_include_directories(kymoknot_qhull PUBLIC
    ${KYMOKNOT_QHULL_DIR}/src/libqhull
)
set_property(TARGET kymoknot_qhull PROPERTY C_STANDARD 99)
set_property(TARGET kymoknot_qhull PROPERTY C_STANDARD_REQUIRED ON)
set_property(TARGET kymoknot_qhull PROPERTY POSITION_INDEPENDENT_CODE ON)
if (UNIX)
    target_link_libraries(kymoknot_qhull PUBLIC m)
endif()

add_library(kymoknot_core STATIC
    ${KYMOKNOT_DIR}/Alexander.c
    ${KYMOKNOT_DIR}/KNT_arc.c
    ${KYMOKNOT_DIR}/KNT_closures.c
    ${KYMOKNOT_DIR}/KNT_qhull.c
    ${KYMOKNOT_DIR}/KNT_simplify.c
    ${KYMOKNOT_DIR}/KNT_identify.c
    ${KYMOKNOT_DIR}/KNT_lib.c
    ${KYMOKNOT_DIR}/lib/my_geom.c
    ${KYMOKNOT_DIR}/lib/my_memory.c
    ${KYMOKNOT_DIR}/lib/messages.c
    ${KYMOKNOT_DIR}/lib/prng.c
    KymoKnotNativeBridge.c
)
target_include_directories(kymoknot_core PRIVATE
    ${KYMOKNOT_DIR}
    ${KYMOKNOT_DIR}/lib
    ${KYMOKNOT_QHULL_DIR}/src/libqhull
)
target_link_libraries(kymoknot_core PUBLIC kymoknot_qhull)
set_property(TARGET kymoknot_core PROPERTY C_STANDARD 99)
set_property(TARGET kymoknot_core PROPERTY C_STANDARD_REQUIRED ON)
set_property(TARGET kymoknot_core PROPERTY POSITION_INDEPENDENT_CODE ON)

add_executable(cablegen
    main.cpp
    CablePhysics.cpp
    CablePhysics.h
    CableGeneration.cpp
    CableGeneration.h
    TanglementAnalyzer.cpp
    TanglementAnalyzer.h
    TopologyAudit.cpp
    TopologyAudit.h
    KymoKnotAnalyzer.cpp
    KymoKnotAnalyzer.h
    KymoKnotNativeBridge.h
)

target_include_directories(cablegen PRIVATE
    ${PROJECT_PATH}
    ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(cablegen PRIVATE
    PositionBasedDynamics
    Simulation
    Utils
    kymoknot_core
)

set_property(TARGET cablegen PROPERTY CXX_STANDARD 17)
set_property(TARGET cablegen PROPERTY CXX_STANDARD_REQUIRED ON)
