CXXFLAGS = -I../include -O2 -Wall

OBJ=numeric.o vector.o matrix.o \
	dense_differential.o sparse_differential.o \
	approximate_taylor_model.o \
	set.o hybrid_set.o \
	point.o box.o zonotope.o \
	grid.o paving.o \
	dynamical_toolbox.o \
	hybrid_automaton.o \
	graphics.o


	#taylor_model.o  

lib: $(OBJ)
	g++ -shared -o libariadne-lite.so $(OBJ) -lplplotcxxd -lcairo -lgtk-x11-2.0

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
	g++ -fpic -I../include -I/usr/include/cairo / -Wall -O2 -c $<

GRAPHICS_CFLAGS  = `pkg-config gtk+-2.0 --cflags`

graphics.o:
	g++ -fpic $(CXXFLAGS) $(GRAPHICS_CFLAGS) -c graphics.cc

include Makefile.dep