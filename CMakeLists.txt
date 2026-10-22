cmake_minimum_required(VERSION 3.12)
project(Bubble VERSION 1.14 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

find_package(GSL 2.0 REQUIRED)
find_package(Threads REQUIRED)

# BUILD_SHARED_LIBS selects static (default) or shared. No CLI entry point is
# compiled into this target, and consumers need only the public header.
add_library(bubble bubble.cc Faddeeva.cc)
add_library(Bubble::bubble ALIAS bubble)
target_compile_definitions(bubble PRIVATE BUBBLE_LIBRARY)
target_compile_features(bubble PUBLIC cxx_std_11)
set_target_properties(bubble PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_include_directories(bubble PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_link_libraries(bubble PUBLIC GSL::gsl GSL::gslcblas Threads::Threads)

set(bubble_package_dir "${CMAKE_INSTALL_LIBDIR}/cmake/Bubble")
install(TARGETS bubble EXPORT BubbleTargets
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(FILES bubble.hh DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
install(EXPORT BubbleTargets FILE BubbleTargets.cmake NAMESPACE Bubble::
  DESTINATION "${bubble_package_dir}")
configure_package_config_file(cmake/BubbleConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/BubbleConfig.cmake"
  INSTALL_DESTINATION "${bubble_package_dir}")
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/BubbleConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}" COMPATIBILITY SameMinorVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/BubbleConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/BubbleConfigVersion.cmake"
  DESTINATION "${bubble_package_dir}")

option(BUBBLE_BUILD_TESTS "Build native Bubble API tests" OFF)
if(BUBBLE_BUILD_TESTS)
  enable_testing()
  add_executable(bubble_api_tests regression/library_tests.cc)
  target_link_libraries(bubble_api_tests PRIVATE Bubble::bubble)
  add_test(NAME bubble_api COMMAND bubble_api_tests "${CMAKE_CURRENT_SOURCE_DIR}")
  foreach(kind IN ITEMS static shared)
    add_test(NAME bubble_installed_${kind}
      COMMAND "${CMAKE_COMMAND}"
        "-DBUBBLE_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}"
        "-DBUBBLE_TEST_DIR=${CMAKE_CURRENT_BINARY_DIR}"
        "-DBUBBLE_GENERATOR=${CMAKE_GENERATOR}"
        "-DBUBBLE_CXX_COMPILER=${CMAKE_CXX_COMPILER}"
        "-DBUBBLE_CTEST_COMMAND=${CMAKE_CTEST_COMMAND}"
        "-DBUBBLE_KIND=${kind}"
        -P "${CMAKE_CURRENT_SOURCE_DIR}/regression/check_installed_package.cmake")
  endforeach()
endif()
