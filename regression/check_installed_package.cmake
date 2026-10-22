cmake_minimum_required(VERSION 3.12)

function(run directory)
  execute_process(COMMAND ${ARGN} WORKING_DIRECTORY "${directory}"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE errors TIMEOUT 120)
  if(NOT "${status}" STREQUAL "0")
    message(FATAL_ERROR "Command failed (${status}): ${ARGN}\n${output}\n${errors}")
  endif()
endfunction()

# Each invocation owns a separate temporary tree, including when CTest is parallel.
string(RANDOM LENGTH 12 ALPHABET 0123456789abcdef nonce)
set(root "${BUBBLE_TEST_DIR}/package-test-${BUBBLE_KIND}-${nonce}")
file(MAKE_DIRECTORY "${root}/producer" "${root}/consumer-build")
set(shared OFF)
set(libdir lib)
set(includedir include)
if(BUBBLE_KIND STREQUAL "shared")
  set(shared ON)
  # Exercise GNUInstallDirs overrides as well as the default layout.
  set(libdir lib64)
  set(includedir include/bubble)
elseif(NOT BUBBLE_KIND STREQUAL "static")
  message(FATAL_ERROR "Unknown package test kind: ${BUBBLE_KIND}")
endif()

run("${root}/producer" "${CMAKE_COMMAND}" "${BUBBLE_SOURCE_DIR}"
  -G "${BUBBLE_GENERATOR}" "-DCMAKE_CXX_COMPILER=${BUBBLE_CXX_COMPILER}"
  -DCMAKE_BUILD_TYPE=Release "-DBUILD_SHARED_LIBS=${shared}" -DBUBBLE_BUILD_TESTS=OFF
  "-DCMAKE_INSTALL_PREFIX=${root}/stage"
  "-DCMAKE_INSTALL_LIBDIR=${libdir}" "-DCMAKE_INSTALL_INCLUDEDIR=${includedir}")
run("${root}/producer" "${CMAKE_COMMAND}" --build . --config Release --target install)

# The consumer has no source-tree include or link paths. Moving the prefix also
# detects hardcoded install locations, while removing the producer rules out
# accidental links to build artifacts.
file(RENAME "${root}/stage" "${root}/installed")
file(REMOVE_RECURSE "${root}/producer")
file(COPY "${BUBBLE_SOURCE_DIR}/regression/installed_consumer/"
  DESTINATION "${root}/consumer-src")
set(package_dir "${root}/installed/${libdir}/cmake/Bubble")
foreach(file IN ITEMS BubbleConfig.cmake BubbleConfigVersion.cmake BubbleTargets.cmake)
  if(NOT EXISTS "${package_dir}/${file}")
    message(FATAL_ERROR "Missing installed package file: ${file}")
  endif()
endforeach()
run("${root}/consumer-build" "${CMAKE_COMMAND}" "${root}/consumer-src"
  -G "${BUBBLE_GENERATOR}" "-DCMAKE_CXX_COMPILER=${BUBBLE_CXX_COMPILER}"
  -DCMAKE_BUILD_TYPE=Release "-DCMAKE_PREFIX_PATH=${root}/installed"
  -DCMAKE_FIND_PACKAGE_NO_PACKAGE_REGISTRY=ON
  -DCMAKE_FIND_PACKAGE_NO_SYSTEM_PACKAGE_REGISTRY=ON)
run("${root}/consumer-build" "${CMAKE_COMMAND}" --build . --config Release)
run("${root}/consumer-build" "${BUBBLE_CTEST_COMMAND}" -C Release --output-on-failure)
message(STATUS "Installed ${BUBBLE_KIND} package verified: ${package_dir}")
