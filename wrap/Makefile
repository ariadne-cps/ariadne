INTERVAL_OBJ=interval_module.o export_numeric.o export_linear_algebra.o 
INTERVAL_LIB=../src/numeric.o ../src/vector.o ../src/matrix.o

DIFFERENTIATION_OBJ=differentiation_module.o export_differentiation.o 
DIFFERENTIATION_LIB=../src/numeric.o ../src/vector.o ../src/matrix.o \
	../src/dense_differential.o ../src/sparse_differential.o 

ARIADNE_OBJ=ariadne_module.o \
	export_numeric.o export_linear_algebra.o export_differentiation.o export_function.o

ARIADNE_LIB=../src/libariadne-lite.so

wrappers: interval.so differentiation.so ariadne.so

dep:  
	@rm -f Makefile.dep
	@for ff in ${OBJ}; \
	do f=`basename $$ff .o`.cc ; \
	g++ -MM -I../include -I/usr/include/python2.5 $$f; echo ""; \
	done >> Makefile.dep

clean:
	rm -f *.o *.so *.dep 

.cc.o:
	g++ -fpic -I../include -I/usr/include/python2.5 -Wall -O2 -c $<

interval.so: $(INTERVAL_OBJ)
	g++ -shared -Wall -O2 -lboost_python $(INTERVAL_OBJ) $(INTERVAL_LIB) -o interval.so

differentiation.so: $(DIFFERENTIATION_OBJ)
	g++ -shared -Wall -O2 -lboost_python $(DIFFERENTIATION_OBJ) $(DIFFERENTIATION_LIB) -o differentiation.so

ariadne.so: ariadne_module.o $(ARIADNE_OBJ)
	g++ -shared -Wall -O2 -lboost_python $(ARIADNE_OBJ) $(ARIADNE_LIB) -o ariadne.so
