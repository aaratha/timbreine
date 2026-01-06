# tests/CMakeLists.txt
cmake_minimum_required(VERSION 3.24)
Include(FetchContent)

set(TEST_EXEC timbreine_tests)

set(TEST_SOURCES
    test_centroid.cpp
    test_flux.cpp
    test_rolloff.cpp
    test_mfcc.cpp
    test_spectral_env.cpp
)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.8.1 # or a later release
)

FetchContent_MakeAvailable(Catch2)

add_executable(${TEST_EXEC} ${TEST_SOURCES})

target_link_libraries(${TEST_EXEC}
    PRIVATE
        Catch2::Catch2WithMain
        timbreine_core
)

target_include_directories(${TEST_EXEC}
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

include(CTest)
include(Catch)

catch_discover_tests(${TEST_EXEC}
    TEST_PREFIX "timbreine::"
)
