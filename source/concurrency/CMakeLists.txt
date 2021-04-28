set(LIBRARY_NAME ariadne-concurrency)

add_library(${LIBRARY_NAME} OBJECT
        buffered_thread.cpp
        thread.cpp
        thread_pool.cpp
        task_manager.cpp
)

foreach(WARN ${LIBRARY_EXCLUSIVE_WARN})
    target_compile_options(${LIBRARY_NAME} PUBLIC -W${WARN})
endforeach()

if(COVERAGE)
    target_link_libraries(${LIBRARY_NAME} PUBLIC coverage_config)
endif()
