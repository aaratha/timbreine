cmake_minimum_required(VERSION 3.24)
project(timbre VERSION 0.1 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Tell CMake where Qt is installed (Homebrew example)
set(CMAKE_PREFIX_PATH "/opt/homebrew/opt/qt6")

# Miniaudio wrapper
set(EXTERNAL_SRC
    "${PROJECT_SOURCE_DIR}/external/miniaudio_wrapper.cpp"
)

# Miniaudio
add_library(miniaudio STATIC external/miniaudio/miniaudio.c)
set_target_properties(miniaudio PROPERTIES LINKER_LANGUAGE C)
target_include_directories(miniaudio PUBLIC external/miniaudio)

# pffft
add_library(pffft STATIC
    external/pffft/pffft.c
    external/pffft/pffft_common.c
)
set_target_properties(pffft PROPERTIES LINKER_LANGUAGE C)
target_include_directories(pffft PUBLIC external/pffft)

# Find required Qt components
find_package(Qt6 REQUIRED COMPONENTS Quick Core Gui Qml)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Standard Qt project setup
qt_standard_project_setup()

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Collect source files
file(GLOB_RECURSE SOURCES
    src/*.cpp
)

# Add executable
qt_add_executable(${PROJECT_NAME}
  ${SOURCES} ${EXTERNAL_SRC} include/ui.hpp
)

qt_add_qml_module(${PROJECT_NAME}
    URI MyApp
    QML_FILES
        qml/MainView.qml
)

# Copy all QML files from source to build folder
file(GLOB_RECURSE QML_SOURCES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} qml/*.qml)

foreach(QML_FILE ${QML_SOURCES})
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${QML_FILE}
                   ${CMAKE_CURRENT_BINARY_DIR}/${QML_FILE}
                   COPYONLY)
endforeach()

# copy example.wav into build directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/example.wav
               ${CMAKE_CURRENT_BINARY_DIR}/example.wav
               COPYONLY)
 

target_compile_definitions(pffft PRIVATE PFFFT_ENABLE_FLOAT)

# Link Qt libraries
target_link_libraries(${PROJECT_NAME}
  PRIVATE
  Qt6::Quick
  Qt6::Core
  Qt6::Gui
  Qt6::Qml
  miniaudio
  pffft
)

# Optional: specify include directories for the target
target_include_directories(${PROJECT_NAME} PRIVATE include external/miniaudio)
