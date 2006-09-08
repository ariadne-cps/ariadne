include ../config.mk

install:
	mkdir -p ${PYTHONPREFIX}/ariadne
	cp -f src/*.py ${PYTHONPREFIX}/ariadne/

check:
	(cd test; make check);


clean:
	rm -rf src/*.pyc test/*.pyc test/*.log test/*.eps


