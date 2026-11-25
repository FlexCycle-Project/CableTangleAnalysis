# MyCableDemo only loads shaders at runtime (geometry is generated
# procedurally), so only the shader and font directories are copied.
add_custom_target(CopyPBDShaders
	${CMAKE_COMMAND} -E copy_directory
	${PROJECT_SOURCE_DIR}/data/shaders
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/resources/shaders
	COMMENT "Copying PBD shaders"
)
set_target_properties(CopyPBDShaders PROPERTIES FOLDER "Data copy")

add_custom_target(CopyImguiFonts
	${CMAKE_COMMAND} -E copy_directory
	${PROJECT_SOURCE_DIR}/extern/imgui/misc/fonts
	${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/resources/fonts
	COMMENT "Copying fonts"
)
set_target_properties(CopyImguiFonts PROPERTIES FOLDER "Data copy")
