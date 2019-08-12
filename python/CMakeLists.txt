add_subdirectory(source)

add_library(pyariadne MODULE
    $<TARGET_OBJECTS:pyariadne-module-obj>
    $<TARGET_OBJECTS:pyariadne-numeric-submodule-obj>
    $<TARGET_OBJECTS:pyariadne-algebra-submodule-obj>
    $<TARGET_OBJECTS:pyariadne-extra-submodule-obj>
)
set_target_properties(pyariadne PROPERTIES PREFIX "" OUTPUT_NAME "ariadne" SUFFIX ".so")
target_link_libraries (pyariadne PUBLIC ariadne ${GCOV_LIBRARIES} PRIVATE pybind11::module)

add_library(pyariadne-numeric MODULE 
    $<TARGET_OBJECTS:pyariadne-numeric-module-obj>
    $<TARGET_OBJECTS:pyariadne-numeric-submodule-obj> 
    $<TARGET_OBJECTS:ariadne-utility> 
    $<TARGET_OBJECTS:ariadne-numeric>
)
set_target_properties(pyariadne-numeric PROPERTIES PREFIX "" OUTPUT_NAME "numeric" SUFFIX ".so")
target_link_libraries (pyariadne-numeric PUBLIC ${GCOV_LIBRARIES} ${MPFR_LIBRARIES} ${GMP_LIBRARIES} PRIVATE pybind11::module)

add_library(pyariadne-algebra MODULE 
    $<TARGET_OBJECTS:pyariadne-algebra-module-obj>
    $<TARGET_OBJECTS:pyariadne-numeric-submodule-obj>    
    $<TARGET_OBJECTS:pyariadne-algebra-submodule-obj>
)
set_target_properties(pyariadne-algebra PROPERTIES PREFIX "" OUTPUT_NAME "algebra" SUFFIX ".so")
target_link_libraries (pyariadne-algebra PUBLIC ${GCOV_LIBRARIES} ariadne PRIVATE pybind11::module)

execute_process(COMMAND ${PYTHON_EXECUTABLE} -m site --user-site OUTPUT_VARIABLE PYARIADNE_INSTALL_DIR)
string(STRIP ${PYARIADNE_INSTALL_DIR} PYARIADNE_INSTALL_DIR)
install(TARGETS pyariadne DESTINATION ${PYARIADNE_INSTALL_DIR})

add_subdirectory(tests)