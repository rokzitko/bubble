cmake_minimum_required(VERSION 3.12)
project(BubbleInstalledConsumer LANGUAGES CXX)

find_package(Bubble 1.14 CONFIG REQUIRED)
get_target_property(bubble_definitions Bubble::bubble INTERFACE_COMPILE_DEFINITIONS)
if("BUBBLE_LIBRARY" IN_LIST bubble_definitions)
  message(FATAL_ERROR "BUBBLE_LIBRARY leaked into the installed interface")
endif()

add_executable(installed_consumer main.cc)
target_link_libraries(installed_consumer PRIVATE Bubble::bubble)
enable_testing()
add_test(NAME installed_consumer COMMAND installed_consumer)
