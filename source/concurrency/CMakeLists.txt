set(LIBRARY_NAME ariadne-concurrency)

add_library(${LIBRARY_NAME} OBJECT
        concurrency_manager.cpp
        task_search_parameter.cpp
        task_search_space.cpp
        task_search_point.cpp
        task_appraisal.cpp
        searchable_configuration.cpp
        configuration_property_path.cpp
        )

cotire(${LIBRARY_NAME} CONFIGURATIONS None)

foreach(WARN ${LIBRARY_EXCLUSIVE_WARN})
    target_compile_options(${LIBRARY_NAME} PUBLIC -W${WARN})
endforeach()

if(COVERAGE)
    target_link_libraries(${LIBRARY_NAME} PUBLIC coverage_config)
endif()
