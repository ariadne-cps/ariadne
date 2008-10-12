INTERVAL_OBJ=../src/numeric.o

DIFFERENTIATION_OBJ=../src/numeric.o ../src/vector.o ../src/matrix.o ../src/dense_differential.o

FUNCTION_OBJ=../src/numeric.o ../src/vector.o ../src/matrix.o ../src/dense_differential.o ../src/approximate_taylor_model.o

GRID_OBJ=../src/numeric.o 

OBJ=interval_module.o differentiation_module.o function_module.o grid_module.o

wrappers: interval.so differentiation.so function.so grid.so

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

interval.so: interval_module.o $(INTERVAL_OBJ)
	g++ -shared -Wall -O2 -lboost_python interval_module.o $(INTERVAL_OBJ) -o interval.so

differentiation.so: differentiation_module.o $(DIFFERENTIATION_OBJ)
	g++ -shared -Wall -O2 -lboost_python differentiation_module.o $(DIFFERENTIATION_OBJ) -o differentiation.so

function.so: function_module.o $(FUNCTION_OBJ)
	g++ -shared -Wall -O2 -lboost_python function_module.o $(FUNCTION_OBJ) -o function.so

grid.so: grid_module.o $(GRID_OBJ)
	g++ -shared -Wall -O2 -lboost_python grid_module.o $(GRID_OBJ) -o grid.so

