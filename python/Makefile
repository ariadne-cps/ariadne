include ../config.mk

install:
	mkdir -p ${PYTHONPREFIX}/ariadne
	cp -f src/*.py ${PYTHONPREFIX}/ariadne/

check:
	(cd test/; for file in *.py; do echo $$file; $$file; done);


clean:
	rm -rf src/*.pyc test/*.pyc test/*.log test/*.eps


