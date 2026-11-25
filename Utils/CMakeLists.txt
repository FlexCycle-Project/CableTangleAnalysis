add_library(Utils
    Hashmap.h
    Timing.h
    Logger.h
    IndexedFaceMesh.cpp
    IndexedFaceMesh.h
    IndexedTetMesh.cpp
    IndexedTetMesh.h
    VolumeIntegration.cpp
    VolumeIntegration.h
)

target_include_directories(Utils PUBLIC
    ${PROJECT_PATH}
    ${EIGEN3_INCLUDE_DIR}
)
