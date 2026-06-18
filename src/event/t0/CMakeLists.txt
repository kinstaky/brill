add_library(brill_dssd_match_event dssd_match_event.cpp)
target_include_directories(brill_dssd_match_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(brill_dssd_match_event PUBLIC ROOT::Tree ROOT::RIO ROOT::Core)

add_library(brill_t0_event t0_event.cpp)
target_include_directories(brill_t0_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(brill_t0_event PUBLIC ROOT::Tree ROOT::RIO ROOT::Core)
