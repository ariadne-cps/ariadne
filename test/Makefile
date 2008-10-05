TESTS=test_differential.o test_sparse_differential.o test_approximate_taylor_model.o 

tests: $(TESTS)
	g++ -L../src test_differential.o -lariadne-lite -o test_differential.exe
	g++ -L../src test_sparse_differential.o -lariadne-lite -o test_sparse_differential.exe
	g++ -L../src test_approximate_taylor_model.o -lariadne-lite -o test_approximate_taylor_model.exe

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

