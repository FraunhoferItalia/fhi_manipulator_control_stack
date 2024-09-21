cmake_minimum_required(VERSION 3.8)
project(robust_controller)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

add_compile_options(-std=c++17 -g)
set(CMAKE_BUILD_TYPE Release)
add_subdirectory(ndof_robust_controller)

set(CONTROLLER_LIBRARY_NAME robust_controller_core)
find_package(catkin REQUIRED COMPONENTS
  fhi_manipulator_control
)
find_package(xbot2 REQUIRED)

catkin_package(
  INCLUDE_DIRS 
    include 
    ndof_robust_controller/include 
    ndof_robust_controller/src/parameterTuning 
    ndof_robust_controller/src/robustController
  LIBRARIES ${CONTROLLER_LIBRARY_NAME}
  CATKIN_DEPENDS 
    fhi_manipulator_control 
  DEPENDS 
    ndof_robust_controller
)

catkin_python_setup()

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
)
# Build
add_library(${CONTROLLER_LIBRARY_NAME} SHARED
  src/controller_core.cpp
)
target_include_directories(${CONTROLLER_LIBRARY_NAME} PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ndof_robust_controller/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ndof_robust_controller/src/parameterTuning>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ndof_robust_controller/src/robustController>
  $<INSTALL_INTERFACE:ndof_robust_controller/include>
  $<INSTALL_INTERFACE:ndof_robust_controller/src/parameterTuning>
  $<INSTALL_INTERFACE:ndof_robust_controller/src/robustController>
)
target_link_libraries(${CONTROLLER_LIBRARY_NAME} PUBLIC ndof_robust_controller controller_core)

if(XBOT2_ENABLE_XENO)
  set_xeno_flags(${CONTROLLER_LIBRARY_NAME})
endif()

# Install
install(
  TARGETS ${CONTROLLER_LIBRARY_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_GLOBAL_BIN_DESTINATION}
)
install(
  DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  FILES_MATCHING PATTERN "*.hpp"
)