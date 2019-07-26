if(${KIRK_FOUND})
    set(EXTERNAL_SRC ${EXTERNAL_SRC} kirk_real.cpp)
endif()

add_library(ariadne-external SHARED ${EXTERNAL_SRC})

target_link_libraries(ariadne-external ariadne)
if(${KIRK_FOUND})
    target_link_libraries(ariadne-external kirk)
endif()
