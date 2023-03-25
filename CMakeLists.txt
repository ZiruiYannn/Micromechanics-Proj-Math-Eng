cmake_minimum_required(VERSION 3.18)
project(PME)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)



### Libraries

add_subdirectory(Lib/eigen-3.4.0)
add_subdirectory(Lib/fftw-3.3.10)
add_subdirectory(Lib/googletest-main)



### Source files

## main executable comes here
# add_executable()



### Test files

enable_testing()

add_executable(IOTest Test/IO_test.cpp)

target_link_libraries(IOTest PUBLIC GTest::gtest_main)
target_link_libraries(IOTest PUBLIC Eigen3::Eigen)
target_include_directories(IOTest PRIVATE Src)

add_executable(algorithmTest Test/algorithmTest.cpp)

target_link_libraries(algorithmTest PUBLIC GTest::gtest_main)
target_link_libraries(algorithmTest PUBLIC Eigen3::Eigen)
target_include_directories(algorithmTest PRIVATE Src)

include(GoogleTest)

gtest_discover_tests(IOTest)
gtest_discover_tests(algorithmTest)













