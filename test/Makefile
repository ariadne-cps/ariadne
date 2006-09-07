include ../config.mk

TESTS = \
  test_array test_binary_word \
  test_lattice \
  test_mpfr test_arithmetic test_function test_interval \
  test_vector test_matrix test_linear_algebra \
  test_point \
  test_rectangle test_parallelotope test_zonotope \
  test_simplex test_polyhedron \
  test_list_set test_grid_set test_partition_tree_set \
  test_epsfstream \
  test_interval_newton \
  test_chainreach \
  test_integration_step 

OBJS=$(addsuffix .o,$(TESTS))

CXXFLAGS += -Wall -I../include/ 
LDLIBS = ../src/libariadne.so -lmpfr -lgmp -lgmpxx -lppl #-lblas

all: ${TESTS}

check: $(TESTS)
	@(for test in $(TESTS); do ./$$test; done);

$(TESTS): %: %.o ../src/libariadne.so
	$(CXX) $@.o $(LDXXFLAGS) $(LDLIBS) -o $@

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
	@for ff in ${TESTS}; do g++ -MM ${CXXFLAGS} $$ff.cc; \
	echo "	-$$""{CXX} $$""{CXXFLAGS} -c $$ff.cc "; \
	echo ""; done >> Makefile

depclean:
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile

### Do Not edit this line
