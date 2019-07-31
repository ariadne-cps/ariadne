if(${KIRK_FOUND})
    set(EXTERNAL_SRC ${EXTERNAL_SRC} kirk_real.cpp)
endif()

if (${EXTERNAL_SRC})
    add_library(ariadne-external SHARED ${EXTERNAL_SRC})
    target_link_libraries(ariadne-external ariadne)
endif()

if(${KIRK_FOUND})
    target_link_libraries(ariadne-external kirk)
endif()
