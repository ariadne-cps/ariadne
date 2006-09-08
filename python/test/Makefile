include ../../config.mk

check:
	@pass=0; fail=0; \
	for test in *.py; do \
	  if python $$test > `basename $$test`.log; then \
	    pass=`expr $$pass + 1`; echo "PASS: $$test"; \
	  else \
	    fail=`expr $$fail + 1`; echo "FAIL: $$test"; \
	  fi; \
	done; \
	if test $$fail -eq 0; then \
	  echo "Passed all of $$pass tests."; \
	else \
	  echo "Failed $$fail of" `expr $$pass + $$fail` "tests."; \
	fi


clean:
	rm -rf src/*.pyc *.pyc *.log *.eps
