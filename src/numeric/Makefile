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

### Do Not edit this line
### Everything below is auto-generated
arithmetic.o: arithmetic.cc ../../include/numeric/arithmetic.h \
  ../../include/numeric/../numeric/numerical_types.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/../numeric/../utility/../base/array.h \
  ../../include/numeric/../numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/numeric/../numeric/../numeric/dyadic.h \
  ../../include/numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/approximation.h \
  ../../include/numeric/../numeric/../numeric/numerical_types.h \
  ../../include/numeric/arithmetic.tpl ../../include/numeric/arithmetic.h \
  ../../include/numeric/numerical_types.h

function.o: function.cc ../../include/numeric/function.h \
  /ufs/collins/include/mpfr.h \
  ../../include/numeric/../numeric/numerical_types.h \
  ../../include/numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/../numeric/../declarations.h \
  ../../include/numeric/../numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/../numeric/../numeric/../utility/../base/array.h \
  ../../include/numeric/../numeric/../numeric/../utility/../base/../declarations.h \
  ../../include/numeric/../numeric/../numeric/dyadic.h \
  ../../include/numeric/../numeric/../utility/stlio.h \
  ../../include/numeric/function.tpl ../../include/numeric/function.h \
  ../../include/numeric/../numeric/integer.h \
  ../../include/numeric/../numeric/arithmetic.h \
  ../../include/numeric/../numeric/../numeric/numerical_types.h \
  ../../include/numeric/../numeric/../numeric/approximation.h \
  ../../include/numeric/../numeric/../numeric/../numeric/numerical_types.h \
  ../../include/numeric/../numeric/approximation.h \
  ../../include/numeric/numerical_types.h ../../include/real_typedef.h \
  ../../include/numeric/numerical_types.h

