add_subdirectory(event)
add_subdirectory(energy_calculator)
add_subdirectory(rebuild)
add_subdirectory(t0)

add_library(config config.cpp)
target_include_directories(config PUBLIC ${PROJECT_SOURCE_DIR} ${ROOT_INCLUDE_DIRS})