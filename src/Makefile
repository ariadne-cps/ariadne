OBJ=numeric.o vector.o matrix.o \
	dense_differential.o sparse_differential.o \
	approximate_taylor_model.o \
	point.o box.o zonotope.o \
	grid.o paving.o \
	dynamical_toolbox.o 

	#taylor_model.o  

lib: $(OBJ)
	g++ -shared -o libariadne-lite.so $(OBJ) -lgmpxx

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
	g++ -fpic -I../include -Wall -O2 -c $<

include Makefile.dep