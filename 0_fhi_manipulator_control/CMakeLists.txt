cmake_minimum_required(VERSION 3.8)
project(fhi_manipulator_control)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

add_compile_options(-std=c++17)
set(CMAKE_BUILD_TYPE Release)

set(CORE_CONTROLLER_LIBRARY_NAME controller_core)

# find dependencies
find_package(catkin REQUIRED COMPONENTS
  fhi_robot_dynamic_library
)
find_package(xbot2 REQUIRED)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES ${CORE_CONTROLLER_LIBRARY_NAME}
  CATKIN_DEPENDS fhi_robot_dynamic_library
  DEPENDS xbot2
)

include_directories(${catkin_INCLUDE_DIRS} include)
# Build
add_library(${CORE_CONTROLLER_LIBRARY_NAME} SHARED
  src/controller_core.cpp
  src/parameters.cpp
)
add_dependencies(${CORE_CONTROLLER_LIBRARY_NAME}
  ${${PROJECT_NAME}_EXPORTED_TARGETS}
  ${catkin_EXPORTED_TARGETS}
)

target_link_libraries(${CORE_CONTROLLER_LIBRARY_NAME} PUBLIC
  ${catkin_LIBRARIES}
)

target_include_directories(${CORE_CONTROLLER_LIBRARY_NAME} PUBLIC
  include
  ${catkin_INCLUDE_DIRS}
)

if(XBOT2_ENABLE_XENO)
  set_xeno_flags(${CORE_CONTROLLER_LIBRARY_NAME})
endif()

# Install
install(
  DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  FILES_MATCHING PATTERN "*.hpp"
)
install(TARGETS ${CORE_CONTROLLER_LIBRARY_NAME}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

install(PROGRAMS 
  scripts/controller_package_generator.py
  DESTINATION lib/${PROJECT_NAME}
)

install(DIRECTORY data
  DESTINATION share/${PROJECT_NAME}
)