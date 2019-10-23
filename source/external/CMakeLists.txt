if(${KIRK_FOUND})
  set(EXTERNAL_SRC ${EXTERNAL_SRC} kirk_real.cpp)
endif()

# ND -- 22/10/2019 -- GLPK and Eigen3 library -- BEGIN
if(GLPK_FOUND)
  add_library(glpk_interface INTERFACE)
  add_custom_target(glpk_interface.headers SOURCES glpk_interface.hpp)
endif()

if(EIGEN3_FOUND)
	add_library(eigen3_interface INTERFACE)
	add_custom_target(eigen3_interface.headers SOURCES eigen3_interface.hpp)
endif()
# ND -- 22/10/2019 -- GLPK and Eigen3 library -- END

if(${EXTERNAL_SRC})
  add_library(ariadne-external SHARED ${EXTERNAL_SRC})
  cotire(ariadne-external CONFIGURATIONS None)
  target_link_libraries(ariadne-external ariadne)
  if(COVERAGE)
    target_link_libraries(ariadne-external PUBLIC coverage_config)
  endif()
endif()

if(${KIRK_FOUND})
  target_link_libraries(ariadne-external kirk)
endif()

