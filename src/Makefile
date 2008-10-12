CXXFLAGS = -fpic -O2 -Wall  -I../include 

OBJ=numeric.o vector.o matrix.o \
	dense_differential.o sparse_differential.o \
	approximate_taylor_model.o \
	set.o hybrid_set.o \
	point.o box.o zonotope.o \
	grid.o paving.o \
	dynamical_toolbox.o \
	hybrid_automaton.o \
	hybrid_evolver.o \
	logging.o graphics.o



lib: $(OBJ)
	g++ -shared -o libariadne-lite.so $(OBJ) `pkg-config gtk+-2.0 --libs`

dep:  
	@rm -f Makefile.dep
	@for ff in ${OBJ}; \
	do f=`basename $$ff .o`.cc ; \
	g++ -MM -I../include $$f; echo ""; \
	done >> Makefile.dep

main: lib main.o
	g++ main.o libariadne-lite.so -o ariadne

clean:
	rm -f *.o *.so *.dep

.cc.o:
	g++ $(CXXFLAGS) -c $<

GRAPHICS_CFLAGS  = `pkg-config gtk+-2.0 --cflags`

graphics.o:
	g++ $(CXXFLAGS) $(GRAPHICS_CFLAGS) -c graphics.cc

plotting: lib plotting.o
	g++ plotting.o libariadne-lite.so `pkg-config gtk+-2.0 --libs` -o plotting

include Makefile.dep