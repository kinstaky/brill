add_subdirectory(ingot)
add_subdirectory(t0)

add_library(brill_particle_event particle_event.cpp)
target_include_directories(brill_particle_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(brill_particle_event PUBLIC ROOT::Tree ROOT::RIO ROOT::Core)

add_library(brill_event INTERFACE)
target_link_libraries(brill_event
	INTERFACE
	brill_tafd_event
	brill_dssd_event
	brill_jbeam_event
	brill_trigger_event
	brill_silicon_event
	# brill_csi_event
	brill_ppac_event
	brill_align_event
	brill_dssd_match_event
	brill_t0_event
	brill_particle_event
)