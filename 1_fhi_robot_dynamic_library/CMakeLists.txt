cmake_minimum_required(VERSION 3.8)
project(fhi_robot_dynamic_library)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

add_compile_options(-std=c++17)
set(CMAKE_BUILD_TYPE Release)

set(LIBRARY_NAME chain)
find_package(catkin REQUIRED COMPONENTS)
find_package(xbot2 REQUIRED)
find_package(yaml-cpp REQUIRED)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${LIBRARY_NAME}
  DEPENDS yaml-cpp
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${YAML_CPP_INCLUDE_DIR}
)

add_library(${LIBRARY_NAME} SHARED
  src/chain.cpp
)

target_link_libraries(${LIBRARY_NAME} ${catkin_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_include_directories(${LIBRARY_NAME}
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Install library
install(
  TARGETS ${LIBRARY_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_GLOBAL_BIN_DESTINATION}
)

if(XBOT2_ENABLE_XENO)
  set_xeno_flags(${LIBRARY_NAME})
endif()

catkin_install_python(PROGRAMS scripts/generate_chain_configuration.py
                      DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION})

add_executable(dyn_functions_test
tests/dyn_functions_test.cpp
)
target_link_libraries(dyn_functions_test ${LIBRARY_NAME})
install(TARGETS 
  dyn_functions_test
  DESTINATION lib/${PROJECT_NAME})

install(
  DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  FILES_MATCHING PATTERN "*.hpp"
)