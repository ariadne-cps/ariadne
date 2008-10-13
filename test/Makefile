TESTS=test_differential.o test_approximate_taylor_model.o \
	test_graphics.o test_serialization.o

tests: $(TESTS)
	g++ -L../src test_differential.o -lariadne-lite -o test_differential.exe
	g++ -L../src test_approximate_taylor_model.o -lariadne-lite -o test_approximate_taylor_model.exe
	g++ -L../src test_serialization.o -lariadne-lite -lboost_serialization -o test_serialization.exe
	g++ -L../src test_graphics.o -lariadne-lite -o test_graphics.exe

dep:  
	@rm -f Makefile.dep
	@for ff in ${TESTS}; \
	do f=`basename $$ff .o`.cc ; \
	g++ -MM -I../include $$f; echo ""; \
	done >> Makefile.dep

clean:
	rm -f *.o *.exe *.dep

.cc.o:
	g++ -fpic -I../include -Wall -O2 -c $<

