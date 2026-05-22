add_executable(match_dssd match_dssd.cpp)
target_link_libraries(
	match_dssd
	PRIVATE
	config
	t0_dssd
)

add_executable(track_t0 track_t0.cpp)
target_link_libraries(
	track_t0
	PRIVATE
	config
	dssd_match_event
	t0_event
	silicon_event
	ROOT::Graf
)

add_executable(rebuild_t0 rebuild_t0.cpp)
target_link_libraries(
	rebuild_t0
	PRIVATE
	config
	t0_event
	particle_event
	energy_calculator
)

add_executable(calibrate_t0 calibrate_t0.cpp)
target_link_libraries(
	calibrate_t0
	PRIVATE
	config
	t0_event
	energy_calculator
	ROOT::Minuit
)

add_subdirectory(estimate)
add_subdirectory(rebuild)
