include ../../config.mk

CXXFLAGS += -I/usr/include/python2.4/
LDXXFLAGS += -L../src
LDLIBS = -lariadne  -lmpfr -lboost_python -lppl ${LIBGMPXX} 

BASE_OBJ=base_module.o \
	export_array.o export_binary_tree.o 

NUMERIC_OBJ=numeric_module.o \
	export_numeric.o export_interval.o export_function.o \

LINEAR_ALGEBRA_OBJ = linear_algebra_module.o \
	export_vector.o export_matrix.o export_linear_program.o

GEOMETRY_OBJ = geometry_module.o \
	export_point.o export_rectangle.o export_parallelotope.o \
	export_zonotope.o export_simplex.o export_polyhedron.o \
	export_lattice_set.o export_list_set.o export_grid_set.o \
	export_partition_tree_set.o 

SYSTEM_OBJ = system_module.o \
	export_lattice_map.o \
	export_map.o export_affine_map.o export_affine_multimap.o \
	export_polynomial_map.o \
	export_vector_field.o export_affine_vector_field.o \
	export_henon_map.o export_lorenz_system.o 

EVALUATION_OBJ = evaluation_module.o \
	export_apply.o \
	export_integrate.o

OUTPUT_OBJ = output_module.o \
	export_postscript_output.o

all: python 

python: base.so numeric.so linear_algebra.so \
	geometry.so system.so evaluation.so \
	output.so

init_file:
	@echo -n '__all__=["base","numeric", '> __init__.py
	@echo '"linear_algebra","geometry","system","evaluation"]' >> __init__.py
	@echo >> __init__.py
	@echo 'from base import *' >> __init__.py
	@echo 'from numeric import *' >> __init__.py
	@echo -n 'Real=' >> __init__.py
	if [ "${REAL_TYPE}" == "RATIONAL" ]; then \
	  echo 'Rational' >> __init__.py; \
	elif [ "${REAL_TYPE}" == "MPFLOAT" ]; then \
	  echo 'MPFloat' >> __init__.py; \
	else \
	  echo "Unknown REAL type"; \
	  false; \
	fi 
	@echo "numeric.Real=Real" >> __init__.py

base.so: ${BASE_OBJ}
	${CXX} -shared -o base.so ${LDXXFLAGS} ${BASE_OBJ} ${LDLIBS} ${RPATHFLAGS}

numeric.so: ${NUMERIC_OBJ}
	${CXX} -shared -o numeric.so ${LDXXFLAGS} ${NUMERIC_OBJ} ${LDLIBS} ${RPATHFLAGS}

linear_algebra.so: ${LINEAR_ALGEBRA_OBJ}
	${CXX} -shared -o linear_algebra.so ${LDXXFLAGS} ${LINEAR_ALGEBRA_OBJ} ${LDLIBS} ${RPATHFLAGS}

geometry.so: ${GEOMETRY_OBJ}
	${CXX} -shared -o geometry.so ${LDXXFLAGS} ${GEOMETRY_OBJ} ${LDLIBS} ${RPATHFLAGS}

system.so: ${SYSTEM_OBJ}
	${CXX} -shared -o system.so ${LDXXFLAGS} ${SYSTEM_OBJ} ${LDLIBS} ${RPATHFLAGS}

evaluation.so: ${EVALUATION_OBJ}
	${CXX} -shared -o evaluation.so ${LDXXFLAGS} ${EVALUATION_OBJ} ${LDLIBS} ${RPATHFLAGS}

output.so: ${OUTPUT_OBJ}
	${CXX} -shared -o output.so ${LDXXFLAGS} ${OUTPUT_OBJ} ${LDLIBS} ${RPATHFLAGS}


install: python init_file  
	mkdir -p ${PYTHONPREFIX}/ariadne
	rm -rf ${PYTHONPREFIX}/ariadne/*.py
	rm -rf ${PYTHONPREFIX}/ariadne/*.pyc
	rm -rf ${PYTHONPREFIX}/ariadne/*.so
	cp -f *.so *.py ${PYTHONPREFIX}/ariadne/

clean:
	rm -f *.o *.so ariadne/*.so

dep:  
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile
	@echo "### Everything below is auto-generated" >> Makefile
	@for ff in ${BASE_OBJ} ${GEOMETRY_OBJ} ${EVALUATION_OBJ} ${SYSTEM_OBJ} ${LINEAR_ALGEBRA_OBJ}; \
	do f=`basename $$ff .o`.cc ; \
	g++ -MM ${CXXFLAGS} -I${PYTHONINCLUDEDIR} $$f; echo ""; \
	done >> Makefile

depclean:
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile

### Do Not edit this line
