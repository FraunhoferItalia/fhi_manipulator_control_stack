cmake_minimum_required(VERSION 3.8)
project(cartesian_impedance_controller)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

add_compile_options(-std=c++17)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

find_package(catkin REQUIRED COMPONENTS
  fhi_manipulator_control
)
find_package(xbot2 REQUIRED)

catkin_package(
  LIBRARIES ${PROJECT_NAME}_core
  CATKIN_DEPENDS 
    fhi_manipulator_control 
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
)
# Build
add_library(${PROJECT_NAME}_core SHARED
  src/controller_core.cpp
)
target_include_directories(${PROJECT_NAME}_core PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(${PROJECT_NAME}_core PUBLIC controller_core)

if(XBOT2_ENABLE_XENO)
  set_xeno_flags(${PROJECT_NAME}_core)
endif()

# Install
install(
  TARGETS ${PROJECT_NAME}_core
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_GLOBAL_BIN_DESTINATION}
)
install(
  DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  FILES_MATCHING PATTERN "*.hpp"
)