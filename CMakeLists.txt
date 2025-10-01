cmake_minimum_required(VERSION 3.21...4.1)
project(navtools)

# --- Version Control ---
set(NAVTOOLS_MAJOR_VERSION 1)
set(NAVTOOLS_MINOR_VERSION 0)
set(NAVTOOLS_PATCH_VERSION 0)
set(NAVTOOLS_VERSION_NUMBER ${NAVTOOLS_MAJOR_VERSION}.${NAVTOOLS_MINOR_VERSION}.${NAVTOOLS_PATCH_VERSION})
message(STATUS "${PROJECT_NAME} version ${NAVTOOLS_VERSION_NUMBER}")

# --- CMake Options ---
option(BUILD_WITH_BLAS "Enable compilation of with BLAS (default is OFF)" OFF)
option(NAVTOOLS_BUILD_TESTS "Enable compilation of navtools tests (default is OFF)" OFF)

# --- Compiler Flags ---
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-ffast-math -DNDEBUG -march=native -03")

# --- Install path configuration ---
include(GNUInstallDirs)
set(INCLUDE_INSTALL_DIR 
    "${CMAKE_INSTALL_INCLUDEDIR}/navtools" 
    CACHE PATH "The directory relative to CMAKE_INSTALL_PREFIX where navtools header files are installed"
)
set(CMAKE_PACKAGE_INSTALL_DIR
    "${CMAKE_INSTALL_DATADIR}/navtools/cmake"
    CACHE PATH "The directory relative to CMAKE_INSTALL_PREFIX where navtoolsConfig.cmake is installed"
)
set(PKGCONFIG_INSTALL_DIR
    "${CMAKE_INSTALL_DATADIR}/pkgconfig"
    CACHE PATH "The directory relative to CMAKE_INSTALL_PREFIX where navtools.pc is installed"
)
foreach(var INCLUDE_INSTALL_DIR CMAKEPACKAGE_INSTALL_DIR PKGCONFIG_INSTALL_DIR)
  # If an absolute path is specified, make it relative to "{CMAKE_INSTALL_PREFIX}".
  if(IS_ABSOLUTE "${${var}}")
    file(RELATIVE_PATH "${var}" "${CMAKE_INSTALL_PREFIX}" "${${var}}")
  endif()
endforeach()

# --- Create the navtools C++ Library ---
find_package(Eigen3 REQUIRED)

add_library(${PROJECT_NAME} INTERFACE)
add_library(navtools::navtools ALIAS ${PROJECT_NAME})
target_include_directories (${PROJECT_NAME} 
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>
)
list(APPEND EXT_LIB "Eigen3::Eigen")
if (BUILD_WITH_BLAS)
    find_package(BLAS REQUIRED)
    list(APPEND EXT_LIB "${BLAS_LIBRARIES}")
    target_compile_definitions(${PROJECT_NAME} INTERFACE EIGEN_USE_BLAS)
endif()
target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen)
set_target_properties(${PROJECT_NAME} PROPERTIES POSITION_INDEPENDENT_CODE ON)

# --- Install rules ---
configure_file(cmake/navtools.pc.in navtools.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/navtools.pc
        DESTINATION ${PKGCONFIG_INSTALL_DIR}
)

# --- Make Library 'Findable' for other CMake Packages ---
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION ${NAVTOOLS_VERSION_NUMBER}
    COMPATIBILITY SameMajorVersion
)
configure_package_config_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/navtoolsConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/navtoolsConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_PACKAGE_INSTALL_DIR}
    NO_SET_AND_CHECK_MACRO 
    NO_CHECK_REQUIRED_COMPONENTS_MACRO 
)

# --- Install package targets ---
install(DIRECTORY ${PROJECT_NAME} 
        DESTINATION ${INCLUDE_INSTALL_DIR} 
        COMPONENT Devel
)
install(TARGETS ${PROJECT_NAME} 
        EXPORT ${PROJECT_NAME}Targets
)
export(TARGETS ${PROJECT_NAME} 
       NAMESPACE ${PROJECT_NAME}:: 
       FILE "${PROJECT_NAME}Targets.cmake"
)
export(PACKAGE ${PROJECT_NAME})
install(EXPORT ${PROJECT_NAME}Targets 
        NAMESPACE ${PROJECT_NAME}:: 
        DESTINATION ${CMAKE_PACKAGE_INSTALL_DIR}
)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_PACKAGE_INSTALL_DIR}
)

# --- Build Executables ---
if (NAVTOOLS_BUILD_TESTS)
    find_package(GTest REQUIRED)
    set(NAVTOOLS_TESTS "")

    file(GLOB testfiles "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp")
    foreach(file ${testfiles})
        get_filename_component(testname ${file} NAME_WE)
        add_executable(${testname} ${file})
        target_link_libraries(${testname} PUBLIC ${PROJECT_NAME} GTest::GTest GTest::Main)
        list(APPEND NAVTOOLS_TESTS ${testname})
    endforeach(file ${testfiles})
    install(TARGETS ${NAVTOOLS_TESTS}
            RUNTIME DESTINATION bin
    )
endif()
