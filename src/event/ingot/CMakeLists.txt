add_library(brill_ppac_event ppac_event.cpp)
target_include_directories(brill_ppac_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(brill_ppac_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)

add_library(brill_tafd_event tafd_event.cpp)
target_include_directories(brill_tafd_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(brill_tafd_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)

add_library(brill_dssd_event dssd_event.cpp)
target_include_directories(brill_dssd_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(brill_dssd_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)

add_library(brill_beam_event beam_event.cpp)
target_include_directories(brill_beam_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(brill_beam_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)

add_library(brill_trigger_event trigger_event.cpp)
target_include_directories(brill_trigger_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(brill_trigger_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)

add_library(brill_silicon_event silicon_event.cpp)
target_include_directories(brill_silicon_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(brill_silicon_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)

# add_library(csi_event csi_event.cpp)
# target_include_directories(csi_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
# target_link_libraries(csi_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)

add_library(brill_align_event align_event.cpp)
target_include_directories(brill_align_event PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})
target_link_libraries(brill_align_event PUBLIC ROOT::Core ROOT::Tree ROOT::RIO)