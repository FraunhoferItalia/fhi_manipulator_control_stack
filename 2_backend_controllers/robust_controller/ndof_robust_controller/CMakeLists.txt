find_package(OpenMP)
find_package(nlohmann_json REQUIRED)

#set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED True)
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

include_directories(include
		src/parameterTuning
		src/robustController)


add_library(ndof_robust_controller
    SHARED
    src/robust_controller.cpp
    src/matlab_interface.cpp
    src/parameterTuning/combineVectorElements.cpp
    src/parameterTuning/fft.cpp
    src/parameterTuning/FFTImplementationCallback.cpp
    src/parameterTuning/linspace.cpp
    src/parameterTuning/parameterTuning.cpp
    src/parameterTuning/parameterTuning_data.cpp
    src/parameterTuning/parameterTuning_initialize.cpp
    src/parameterTuning/parameterTuning_terminate.cpp
    src/robustController/minOrMax.cpp
    src/robustController/prodMatMatInt.cpp
    src/robustController/recursiveNEurdf.cpp
    src/robustController/recursiveNEurdfIA.cpp
    src/robustController/robustControllerURDF.cpp
    src/robustController/robustControllerURDF_data.cpp
    src/robustController/robustControllerURDF_initialize.cpp
    src/robustController/robustControllerURDF_terminate.cpp
)

target_link_libraries(ndof_robust_controller nlohmann_json::nlohmann_json OpenMP::OpenMP_CXX)

install(TARGETS ndof_robust_controller
  EXPORT ndof_robust_controller
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

install(EXPORT ndof_robust_controller FILE ndof_robust_controller-targets.cmake DESTINATION share/cmake/foo)

install(DIRECTORY include/
  DESTINATION include/
)

add_executable(ndof_robust_controller_test
    src/controller_test.cpp
)

target_link_libraries(ndof_robust_controller_test ndof_robust_controller)
