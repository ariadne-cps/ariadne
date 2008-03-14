#!/bin/bash
for PY in *.py ; do
  echo "Running " $PY "..." ;
  python $PY ;
done 
