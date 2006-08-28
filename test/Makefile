include ../config.mk

TESTS = \
  test_array test_binary_word \
  test_lattice \
  test_arithmetic test_function test_interval \
  test_vector test_matrix test_linear_algebra \
  test_point \
  test_rectangle test_parallelotope test_zonotope \
  test_simplex test_polyhedron \
  test_list_set test_grid_set test_partition_tree_set \
  test_epsfstream \
  test_interval_newton \
  test_chainreach \
  test_integration_step 

CXXFLAGS += -Wall -I../include/ 
#LDXXFLAGS += .-L../src/ -lariadne -lmpfr -lgmp -lgmpxx -lppl #-lblas
LDXXFLAGS += ../src/libariadne.so -lmpfr -lgmp -lgmpxx -lppl #-lblas

all: ${TESTS}

.cc.o:
	${CXX} ${INCLUDE} ${CXXFLAGS} -c $<

clean:
	-rm -rf core *.log *.eps *.o ${TESTS} 

dep:  
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile
	@echo "### Everything below is auto-generated" >> Makefile
	@echo "" >> Makefile 
	@echo "all: $$""{TESTS}" >> Makefile
	@for ff in ${TESTS}; do echo "	-./$$ff"; done >> Makefile
	@echo "" >> Makefile 
	@for ff in ${TESTS}; do g++ -MM ${CXXFLAGS} $$ff.cc; \
	echo "	-$$""{CXX} $$""{CXXFLAGS} -c $$ff.cc "; \
	echo " "; \
	echo "$$ff: $$ff.o"; \
	echo "	-$$""{CXX} $$ff.o $$""{LDXXFLAGS} -o $$ff "; \
	echo ""; done >> Makefile

depclean:
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile

### Do Not edit this line
