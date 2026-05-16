add_executable(match_dssd match_dssd.cpp)
target_link_libraries(
	match_dssd
	PRIVATE
	estimate_support
)

add_subdirectory(estimate)
