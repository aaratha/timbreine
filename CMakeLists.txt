cmake_minimum_required(VERSION 3.24)
project(timbreine VERSION 0.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Qt (Homebrew example)
set(CMAKE_PREFIX_PATH "/opt/homebrew/opt/qt6")

# ------------------------------------------------------------
# FetchContent
# ------------------------------------------------------------
include(FetchContent)

# ----------------------
# miniaudio
# ----------------------
FetchContent_Declare(
    miniaudio_src
    GIT_REPOSITORY https://github.com/mackron/miniaudio.git
    GIT_TAG master
)

FetchContent_MakeAvailable(miniaudio_src)

# ----------------------
# pffft
# ----------------------
FetchContent_Declare(
    pffft_src
    GIT_REPOSITORY https://github.com/marton78/pffft.git
    GIT_TAG master
)

FetchContent_Populate(pffft_src)

add_library(pffft STATIC
    ${pffft_src_SOURCE_DIR}/pffft.c
    ${pffft_src_SOURCE_DIR}/pffft_common.c
)

target_include_directories(pffft PUBLIC
    ${pffft_src_SOURCE_DIR}
)

target_compile_definitions(pffft PRIVATE PFFFT_ENABLE_FLOAT)

# ----------------------
# umappp
# ----------------------
FetchContent_Declare(
  umappp 
  GIT_REPOSITORY https://github.com/libscran/umappp
  GIT_TAG master # replace with a pinned release
)

FetchContent_MakeAvailable(umappp)

# ------------------------------------------------------------
# Qt
# ------------------------------------------------------------
find_package(Qt6 REQUIRED COMPONENTS Quick Core Gui Qml)

qt_standard_project_setup()

# ------------------------------------------------------------
# Sources
# ------------------------------------------------------------
set(CORE_SOURCES
    src/analysis.cpp
    src/audio.cpp
    src/dsp.cpp
)

set(UI_SOURCES
    src/ui.cpp
)

set(APP_SOURCES
    src/main.cpp
)

add_library(timbreine_core STATIC
    ${CORE_SOURCES}
)

target_include_directories(timbreine_core PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${miniaudio_src_SOURCE_DIR}
)

target_link_libraries(timbreine_core PUBLIC
    miniaudio
    pffft
    libscran::umappp
)

add_library(timbreine_ui STATIC
    ${UI_SOURCES}
    include/ui.hpp
)

target_include_directories(timbreine_ui PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(timbreine_ui PUBLIC
    timbreine_core
    Qt6::Core
)

qt_add_executable(${PROJECT_NAME}
    ${APP_SOURCES}
)

qt_add_qml_module(${PROJECT_NAME}
    URI MyApp
    QML_FILES
        qml/MainView.qml
)

# ------------------------------------------------------------
# QML copy (optional but fine)
# ------------------------------------------------------------
file(GLOB_RECURSE QML_SOURCES RELATIVE
    ${CMAKE_CURRENT_SOURCE_DIR}
    qml/*.qml
)

foreach(QML_FILE ${QML_SOURCES})
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/${QML_FILE}
        ${CMAKE_CURRENT_BINARY_DIR}/${QML_FILE}
        COPYONLY
    )
endforeach()

# ------------------------------------------------------------
# Includes + Linking
# ------------------------------------------------------------
target_link_libraries(${PROJECT_NAME}
    PRIVATE
    Qt6::Quick
    Qt6::Core
    Qt6::Gui
    Qt6::Qml
    timbreine_core
    timbreine_ui
)

# ------------------------------------------------------------
# Tests
# ------------------------------------------------------------
include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
