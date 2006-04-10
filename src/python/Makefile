include ../../config.mk

CXXFLAGS += -O2 -fpic -Wall -I/usr/include/python2.4/
LDXXFLAGS += -L../src
LDLIBS = -lariadne -lboost_python ${LIBGMPXX} -lppl

BASE_OBJ=base_module.o \
	export_numeric.o export_interval.o export_function.o \
	export_array.o export_binary_tree.o 

LINEAR_ALGEBRA_OBJ = linear_algebra_module.o \
	export_vector.o export_matrix.o export_linear_program.o

GEOMETRY_OBJ = geometry_module.o \
	export_point.o export_rectangle.o export_parallelotope.o \
	export_zonotope.o export_simplex.o export_polyhedron.o \
	export_lattice_set.o \
	export_list_set.o export_grid_set.o export_partition_tree_set.o \
	export_postscript_output.o \

EVALUATION_OBJ = evaluation_module.o \
	export_map.o export_affine_map.o export_polynomial_map.o \
	export_vector_field.o export_affine_vector_field.o \
	export_apply.o \
	export_integrate.o \
	export_henon_map.o export_lorenz_system.o

all: python

python: base.so linear_algebra.so geometry.so evaluation.so 

base.so: ${BASE_OBJ}
	${CXX} -shared -o base.so ${LDXXFLAGS} ${BASE_OBJ} ${LDLIBS} ${RPATHFLAGS}

linear_algebra.so: ${LINEAR_ALGEBRA_OBJ}
	${CXX} -shared -o linear_algebra.so ${LDXXFLAGS} ${LINEAR_ALGEBRA_OBJ} ${LDLIBS} ${RPATHFLAGS}

geometry.so: ${GEOMETRY_OBJ}
	${CXX} -shared -o geometry.so ${LDXXFLAGS} ${GEOMETRY_OBJ} ${LDLIBS} ${RPATHFLAGS}

evaluation.so: ${EVALUATION_OBJ}
	${CXX} -shared -o evaluation.so ${LDXXFLAGS} ${EVALUATION_OBJ} ${LDLIBS} ${RPATHFLAGS}


install: python
	mkdir -p ${PYTHONINSTALLDIR}/ariadne
	rm -rf ${PYTHONINSTALLDIR}/ariadne/*.pyc
	cp -f *.so *.py ${PYTHONINSTALLDIR}/ariadne/

clean:
	rm -f *.o *.so ariadne/*.so

dep:  
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile
	@echo "### Everything below is auto-generated" >> Makefile
	@for ff in ${BASE_OBJ} ${GEOMETRY_OBJ} ${EVALUATION_OBJ} ${LINEAR_ALGEBRA_OBJ}; \
	do f=`basename $$ff .o`.cc ; \
	g++ -MM ${CXXFLAGS} -I${PYTHONINCLUDEDIR} $$f; echo ""; \
	done >> Makefile

### Do Not edit this line
### Everything below is auto-generated
