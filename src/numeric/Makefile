include ../../config.mk

OBJ = arithmetic.o function.o

all: ${OBJ}

.cc.o: 
	${CXX} ${INCLUDE} ${CXXFLAGS} -c $<

clean:                                                                          
	rm -f *.bak *.o *~ *% core libariadne.*

dep:  
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile
	@echo "### Everything below is auto-generated" >> Makefile
	@for ff in ${OBJ}; \
	do f=`basename $$ff .o`.cc ; \
	g++ -MM ${CXXFLAGS} $$f; echo ""; \
	done >> Makefile

depclean:
	@sed -e '/^### Do Not edit this line$$/,$$d' Makefile \
		> MMM.$$$$ && mv MMM.$$$$ Makefile
	@echo "### Do Not edit this line" >> Makefile

### Do Not edit this line
### Everything below is auto-generated
arithmetic.o: arithmetic.cc ../../include/numeric/approximation.h \
  ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/approximation.tpl \
  ../../include/numeric/arithmetic.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/arithmetic.tpl ../../include/numeric/arithmetic.h \
  ../../include/numeric/float64.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/mpfloat.h \
  ../../include/numeric/../numeric/rational.h \
  ../../include/numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/dyadic.h ../../include/numeric/../base/dyadic.h \
  ../../include/numeric/rational.h

function.o: function.cc ../../include/numeric/function.h \
  ../../include/numeric/integer.h ../../include/numeric/../declarations.h \
  ../../include/numeric/../numeric/numerical_traits.h \
  ../../include/numeric/function.tpl ../../include/numeric/function.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/numerical_types.h \
  ../../include/numeric/../numeric/mpfloat.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/function.h \
  ../../include/numeric/../numeric/../numeric/rational.h \
  ../../include/numeric/../numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/../numeric/../numeric/numerical_traits.h \
  ../../include/numeric/../numeric/../numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../numeric/../numeric/function.h \
  ../../include/numeric/../numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/approximation.h \
  ../../include/numeric/numerical_types.h ../../include/numeric/float64.h \
  ../../include/numeric/../numeric/function.h \
  ../../include/numeric/mpfloat.h ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h

