add_library(PositionBasedDynamics
    ${PROJECT_PATH}/Common/Common.h

    DirectPositionBasedSolverForStiffRodsInterface.h
    MathFunctions.cpp
    MathFunctions.h
    PositionBasedDynamics.cpp
    PositionBasedDynamics.h
    PositionBasedElasticRods.cpp
    PositionBasedElasticRods.h
    PositionBasedFluids.cpp
    PositionBasedFluids.h
    PositionBasedRigidBodyDynamics.cpp
    PositionBasedRigidBodyDynamics.h
    PositionBasedGenericConstraints.h
    SPHKernels.cpp
    SPHKernels.h
    TimeIntegration.cpp
    TimeIntegration.h
    XPBD.cpp
    XPBD.h
)

target_include_directories(PositionBasedDynamics PUBLIC
    ${PROJECT_PATH}
    ${EIGEN3_INCLUDE_DIR}
)
