set(LIBRARY_NAME ariadne-concurrency)

add_library(${LIBRARY_NAME} OBJECT
        buffered_smart_thread.cpp
        thread_pool_worker.cpp
        smart_thread_pool.cpp
        concurrency_manager.cpp
)

foreach(WARN ${LIBRARY_EXCLUSIVE_WARN})
    target_compile_options(${LIBRARY_NAME} PUBLIC -W${WARN})
endforeach()

if(COVERAGE)
    target_link_libraries(${LIBRARY_NAME} PUBLIC coverage_config)
endif()
