add_executable(rebuild_2alpha rebuild_2alpha.cpp)
target_link_libraries(
	rebuild_2alpha
	PRIVATE
	config
	particle_event
	energy_calculator
	rebuild_physics
	ROOT::Hist
)

add_executable(rebuild_3alpha rebuild_3alpha.cpp)
target_link_libraries(
	rebuild_3alpha
	PRIVATE
	config
	particle_event
	energy_calculator
	rebuild_physics
	ROOT::Hist
)
