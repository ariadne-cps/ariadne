cmake_minimum_required(VERSION 3.15)
message(STATUS "CMake version ${CMAKE_VERSION}")

project(Ariadne VERSION 2.5.3)

include(CMakePackageConfigHelpers)

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message(WARNING "Intel compiler not supported.")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    message(WARNING "Microsoft compiler not supported.")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.0)
        message(FATAL_ERROR "GCC version must be at least 12!")
    endif()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if ((NOT APPLE) AND (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 16.0))
        message(FATAL_ERROR "Clang version must be at least 16!")
    endif()
else()
    message(WARNING "You are using an unidentified compiler! Compilation has only been tested with Clang and GCC.")
endif()

set(ANY_TARGET_WARN all extra pedantic sign-conversion cast-qual disabled-optimization
                  init-self missing-include-dirs sign-promo switch-default undef redundant-decls)
set(ANY_TARGET_DONT_WARN strict-aliasing unused-parameter)
set(LIBRARY_EXCLUSIVE_WARN shadow)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    "${CMAKE_SOURCE_DIR}/cmake/"
    "${CMAKE_SOURCE_DIR}/submodules/configuration/cmake")
set(CMAKE_MACOSX_RPATH 1)
set(BUILD_SHARED_LIBS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors" OFF)

foreach(WARN ${ANY_TARGET_WARN})
    add_compile_options(-W${WARN})
endforeach()
foreach(WARN ${ANY_TARGET_DONT_WARN})
    add_compile_options(-Wno-${WARN})
endforeach()
if(WARNINGS_AS_ERRORS)
    add_compile_options(-Werror)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_compile_options(-Wno-non-template-friend)
endif()
if((CMAKE_CXX_COMPILER_ID MATCHES "GNU") OR (NOT APPLE))
    add_compile_options(-frounding-math)
endif()
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    add_compile_options(-fdiagnostics-color=always)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-fcolor-diagnostics)
endif ()

option(COVERAGE "Enable coverage reporting" OFF)

if(COVERAGE)
    if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang")
        include(LLVMCodeCoverage)
        append_llvm_coverage_compiler_flags()
    elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        include(GCCCodeCoverage)
        append_gcc_coverage_compiler_flags()
    else()
        message(FATAL_ERROR "Coverage is supported only with GCC or macOS AppleClang/Clang.")
    endif()
endif()

enable_testing()
include(CTest)

include(UseDoxygen OPTIONAL)

include(CheckTypeSize)
include(CheckIncludeFile)
include(CheckIncludeFiles)
include(CheckIncludeFileCXX)
include(CheckLibraryExists)

check_type_size("char" SIZEOF_CHAR)
check_type_size("short int" SIZEOF_SHORT_INT)
check_type_size("int" SIZEOF_INT)
check_type_size("long int" SIZEOF_LONG_INT)
check_type_size("long long int" SIZEOF_LONG_LONG_INT)
check_type_size("float" SIZEOF_FLOAT)
check_type_size("double" SIZEOF_DOUBLE)

check_include_file("dlfcn.h" HAVE_DLFCN_H)
check_include_file("execinfo.h" HAVE_EXECINFO_H)
check_include_file_cxx("cxxabi.h" HAVE_CXXABI_H)
if((${HAVE_DLFCN_H}) AND (${HAVE_EXECINFO_H}) AND (${HAVE_CXXABI_H}))
    set(ARIADNE_ENABLE_STACK_TRACE ON)
endif()

option(PACKAGED_SUBMODULES "Notify that the submodules are already available in their subdirectories; for release purposes" OFF)

if(NOT PACKAGED_SUBMODULES)
find_package(Git)
if(GIT_FOUND)
    if (NOT EXISTS ${PROJECT_SOURCE_DIR}/.git) # Manages the case when an archive is used
        execute_process(COMMAND git init WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
        execute_process(COMMAND rm -Rf ${PROJECT_SOURCE_DIR}/python/pybind11 OUTPUT_QUIET ERROR_QUIET)
        execute_process(COMMAND git submodule add -f https://github.com/pybind/pybind11 python/pybind11 WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
        execute_process(COMMAND rm -Rf ${PROJECT_SOURCE_DIR}/submodules/helper OUTPUT_QUIET ERROR_QUIET)
        execute_process(COMMAND git submodule add -f https://github.com/ariadne-cps/helper submodules/helper WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
        execute_process(COMMAND rm -Rf ${PROJECT_SOURCE_DIR}/submodules/logging OUTPUT_QUIET ERROR_QUIET)
        execute_process(COMMAND git submodule add -f https://github.com/ariadne-cps/logging submodules/logging WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
        execute_process(COMMAND rm -Rf ${PROJECT_SOURCE_DIR}/submodules/threading OUTPUT_QUIET ERROR_QUIET)
        execute_process(COMMAND git submodule add -f https://github.com/ariadne-cps/threading submodules/threading WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
        execute_process(COMMAND rm -Rf ${PROJECT_SOURCE_DIR}/submodules/configuration OUTPUT_QUIET ERROR_QUIET)
        execute_process(COMMAND git submodule add -f https://github.com/ariadne-cps/configuration submodules/configuration WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
    endif()
    message(STATUS "Getting all Ariadne submodules if not downloaded already...")
    execute_process(COMMAND git submodule update --init --recursive WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
    if(NOT EXISTS ${PROJECT_SOURCE_DIR}/python/pybind11/CMakeLists.txt)
        message(ERROR "pybind11 download failure")
    else()
        message(STATUS "pybind11 submodule found")
    endif()
    if(NOT EXISTS ${PROJECT_SOURCE_DIR}/submodules/helper/CMakeLists.txt)
        message(ERROR "Helper download failure")
    else()
        message(STATUS "Helper submodule found")
    endif()
    if(NOT EXISTS ${PROJECT_SOURCE_DIR}/submodules/logging/CMakeLists.txt)
        message(ERROR "Logging download failure")
    else()
        message(STATUS "Logging submodule found")
    endif()
    if(NOT EXISTS ${PROJECT_SOURCE_DIR}/submodules/threading/CMakeLists.txt)
        message(ERROR "Threading download failure")
    else()
        message(STATUS "Threading submodule found")
    endif()
    if(NOT EXISTS ${PROJECT_SOURCE_DIR}/submodules/configuration/cmake)
        message(ERROR "Configuration download failure")
    else()
        message(STATUS "Configuration submodule found")
    endif()
else()
    message(FATAL_ERROR "Git not found, no submodules (including pybind11 for Python bindings) could be downloaded")
endif()
endif()


if(EXISTS "${PROJECT_SOURCE_DIR}/.git")
    function(ariadne_gitlink_commit REPOSITORY SUBMODULE_PATH OUTPUT_VARIABLE)
        execute_process(
            COMMAND "${GIT_EXECUTABLE}" -C "${REPOSITORY}" rev-parse "HEAD:${SUBMODULE_PATH}"
            OUTPUT_VARIABLE COMMIT
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_VARIABLE ERROR_TEXT
            RESULT_VARIABLE RESULT
        )
        if(NOT RESULT EQUAL 0 OR COMMIT STREQUAL "")
            message(FATAL_ERROR
                "Cannot determine gitlink for ${REPOSITORY}/${SUBMODULE_PATH}: ${ERROR_TEXT}")
        endif()
        set(${OUTPUT_VARIABLE} "${COMMIT}" PARENT_SCOPE)
    endfunction()

    ariadne_gitlink_commit("${PROJECT_SOURCE_DIR}" "submodules/helper" ARIADNE_HELPER_COMMIT)
    ariadne_gitlink_commit("${PROJECT_SOURCE_DIR}" "submodules/logging" ARIADNE_LOGGING_COMMIT)
    ariadne_gitlink_commit("${PROJECT_SOURCE_DIR}" "submodules/configuration" ARIADNE_CONFIGURATION_COMMIT)

    ariadne_gitlink_commit("${PROJECT_SOURCE_DIR}/submodules/threading" "submodules/helper" THREADING_HELPER_COMMIT)
    ariadne_gitlink_commit("${PROJECT_SOURCE_DIR}/submodules/threading" "submodules/logging" THREADING_LOGGING_COMMIT)
    ariadne_gitlink_commit("${PROJECT_SOURCE_DIR}/submodules/threading" "submodules/configuration" THREADING_CONFIGURATION_COMMIT)

    ariadne_gitlink_commit("${PROJECT_SOURCE_DIR}/submodules/helper" "submodules/configuration" HELPER_CONFIGURATION_COMMIT)
    ariadne_gitlink_commit("${PROJECT_SOURCE_DIR}/submodules/logging" "submodules/configuration" LOGGING_CONFIGURATION_COMMIT)

    if(NOT ARIADNE_HELPER_COMMIT STREQUAL THREADING_HELPER_COMMIT)
        message(FATAL_ERROR
            "Helper dependency mismatch:\n"
            "  Ariadne:   ${ARIADNE_HELPER_COMMIT}\n"
            "  Threading: ${THREADING_HELPER_COMMIT}")
    endif()
    if(NOT ARIADNE_LOGGING_COMMIT STREQUAL THREADING_LOGGING_COMMIT)
        message(FATAL_ERROR
            "Logging dependency mismatch:\n"
            "  Ariadne:   ${ARIADNE_LOGGING_COMMIT}\n"
            "  Threading: ${THREADING_LOGGING_COMMIT}")
    endif()
    if(NOT ARIADNE_CONFIGURATION_COMMIT STREQUAL THREADING_CONFIGURATION_COMMIT OR
       NOT ARIADNE_CONFIGURATION_COMMIT STREQUAL HELPER_CONFIGURATION_COMMIT OR
       NOT ARIADNE_CONFIGURATION_COMMIT STREQUAL LOGGING_CONFIGURATION_COMMIT)
        message(FATAL_ERROR
            "Configuration dependency mismatch:\n"
            "  Ariadne:   ${ARIADNE_CONFIGURATION_COMMIT}\n"
            "  Threading: ${THREADING_CONFIGURATION_COMMIT}\n"
            "  Helper:    ${HELPER_CONFIGURATION_COMMIT}\n"
            "  Logging:   ${LOGGING_CONFIGURATION_COMMIT}")
    endif()
    message(STATUS "Verified Helper commit: ${ARIADNE_HELPER_COMMIT}")
    message(STATUS "Verified Logging commit: ${ARIADNE_LOGGING_COMMIT}")
    message(STATUS "Verified shared Configuration commit: ${ARIADNE_CONFIGURATION_COMMIT}")
endif()

find_package(Threads)
if(NOT Threads_FOUND)
    message(FATAL_ERROR "Unable to find Threads library")
endif()

find_package(GMP)
if(GMP_FOUND)
    include_directories(SYSTEM ${GMP_INCLUDES})
else()
    message(FATAL_ERROR "Unable to find GMP")
endif()

find_package(MPFR)
if(MPFR_FOUND)
    include_directories(SYSTEM ${MPFR_INCLUDES})
else()
    message(FATAL_ERROR "Unable to find MPFR")
endif()

find_package(Cairo)
if(CAIRO_FOUND)
    set(HAVE_CAIRO_H ON)
    include_directories(SYSTEM ${CAIRO_INCLUDE_DIRS})
endif()

find_package(Gnuplot 5.2.8)
if(GNUPLOT_FOUND)
    set(HAVE_GNUPLOT_H ON)
endif()

if(NOT CAIRO_FOUND AND NOT GNUPLOT_FOUND)
    message(STATUS "No Cairo or Gnuplot backend support found: figures can't be written to file")
endif()

add_custom_target(everything)
add_custom_target(tests)

configure_file(${PROJECT_SOURCE_DIR}/cmake/config.hpp.cmake ${PROJECT_SOURCE_DIR}/source/config.hpp)

include_directories(
    ${PROJECT_SOURCE_DIR}/source
    ${PROJECT_SOURCE_DIR}/submodules/helper/include
    ${PROJECT_SOURCE_DIR}/submodules/logging/include
    ${PROJECT_SOURCE_DIR}/submodules/threading/include
)

add_subdirectory(submodules)

add_subdirectory(source)
add_subdirectory(python)

add_subdirectory(tests EXCLUDE_FROM_ALL)
add_subdirectory(examples EXCLUDE_FROM_ALL)
add_subdirectory(tutorials EXCLUDE_FROM_ALL)
add_subdirectory(experimental EXCLUDE_FROM_ALL)

add_library(ariadne-core SHARED
    $<TARGET_OBJECTS:ariadne-utility>
    $<TARGET_OBJECTS:ariadne-foundations>
    $<TARGET_OBJECTS:ariadne-numeric>
    $<TARGET_OBJECTS:ariadne-algebra>
    $<TARGET_OBJECTS:THREADING_SRC>
    $<TARGET_OBJECTS:LOGGING_SRC>
    $<TARGET_OBJECTS:HELPER_SRC>
    source/geometry/interval.cpp
)

add_library(ariadne-kernel SHARED
    $<TARGET_OBJECTS:ariadne-utility>
    $<TARGET_OBJECTS:ariadne-foundations>
    $<TARGET_OBJECTS:ariadne-numeric>
    $<TARGET_OBJECTS:ariadne-algebra>
    $<TARGET_OBJECTS:ariadne-function>
    $<TARGET_OBJECTS:ariadne-geometry>
    $<TARGET_OBJECTS:ariadne-solvers>
    $<TARGET_OBJECTS:ariadne-io>
    $<TARGET_OBJECTS:ariadne-symbolic>
    $<TARGET_OBJECTS:THREADING_SRC>
    $<TARGET_OBJECTS:LOGGING_SRC>
    $<TARGET_OBJECTS:HELPER_SRC>
)

add_library(ariadne SHARED
    $<TARGET_OBJECTS:ariadne-utility>
    $<TARGET_OBJECTS:ariadne-foundations>
    $<TARGET_OBJECTS:ariadne-numeric>
    $<TARGET_OBJECTS:ariadne-algebra>
    $<TARGET_OBJECTS:ariadne-function>
    $<TARGET_OBJECTS:ariadne-geometry>
    $<TARGET_OBJECTS:ariadne-solvers>
    $<TARGET_OBJECTS:ariadne-io>
    $<TARGET_OBJECTS:ariadne-symbolic>
    $<TARGET_OBJECTS:ariadne-dynamics>
    $<TARGET_OBJECTS:ariadne-hybrid>
    $<TARGET_OBJECTS:THREADING_SRC>
    $<TARGET_OBJECTS:LOGGING_SRC>
    $<TARGET_OBJECTS:HELPER_SRC>
)

target_include_directories(ariadne PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>
    $<INSTALL_INTERFACE:include/ariadne>
)
target_link_libraries(ariadne PUBLIC GMP::GMP MPFR::MPFR)
target_link_libraries(ariadne PRIVATE dl ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(ariadne-kernel dl ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(ariadne-core dl ${GMP_LIBRARIES} ${MPFR_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

if(COVERAGE)
    if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang")
        setup_target_for_coverage_llvm(
            NAME coverage
            TARGET ariadne
            DEPENDENCIES tests
            EXCLUDE_REGEX ".*/(tests|python|submodules)/.*"
        )
    else()
        setup_target_for_coverage_gcc(
            NAME coverage
            DEPENDENCIES tests
        )
    endif()
endif()

if(CAIRO_FOUND)
    target_link_libraries(ariadne PUBLIC Cairo::Cairo)
    target_link_libraries(ariadne-kernel ${CAIRO_LIBRARIES})
endif()

configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/AriadneConfig.cmake.in
    ${PROJECT_BINARY_DIR}/AriadneConfig.cmake
    INSTALL_DESTINATION lib/cmake/Ariadne
)
write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/AriadneConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(TARGETS ariadne EXPORT AriadneTargets DESTINATION lib)
install(EXPORT AriadneTargets
    FILE AriadneTargets.cmake
    NAMESPACE Ariadne::
    DESTINATION lib/cmake/Ariadne
)
install(FILES
    ${PROJECT_BINARY_DIR}/AriadneConfig.cmake
    ${PROJECT_BINARY_DIR}/AriadneConfigVersion.cmake
    ${PROJECT_SOURCE_DIR}/cmake/FindGMP.cmake
    ${PROJECT_SOURCE_DIR}/cmake/FindMPFR.cmake
    ${PROJECT_SOURCE_DIR}/cmake/FindCairo.cmake
    DESTINATION lib/cmake/Ariadne
)
install(DIRECTORY source/ DESTINATION include/ariadne FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY submodules/helper/include/helper DESTINATION include FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY submodules/logging/include/logging DESTINATION include FILES_MATCHING PATTERN "*.hpp")
install(DIRECTORY submodules/threading/include/threading DESTINATION include FILES_MATCHING PATTERN "*.hpp")
