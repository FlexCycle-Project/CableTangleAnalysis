add_library(Simulation
    AABB.h
    CollisionDetection.cpp
    CollisionDetection.h
    Constraints.cpp
    Constraints.h
    IDFactory.cpp
    IDFactory.h
    LineModel.cpp
    LineModel.h
    NeighborhoodSearchSpatialHashing.cpp
    NeighborhoodSearchSpatialHashing.h
    ParticleData.h
    RigidBody.h
    RigidBodyGeometry.cpp
    RigidBodyGeometry.h
    Simulation.cpp
    Simulation.h
    SimulationModel.cpp
    SimulationModel.h
    TetModel.cpp
    TetModel.h
    TimeManager.cpp
    TimeManager.h
    TimeStep.cpp
    TimeStep.h
    TimeStepController.cpp
    TimeStepController.h
    TriangleModel.cpp
    TriangleModel.h

    BoundingSphere.h
    BoundingSphereHierarchy.cpp
    BoundingSphereHierarchy.h
    kdTree.h
    kdTree.inl
)

target_include_directories(Simulation PUBLIC
    ${PROJECT_PATH}
    ${EIGEN3_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/extern/GenericParameters/include
)

target_link_libraries(Simulation PUBLIC
    PositionBasedDynamics
    Utils
)
