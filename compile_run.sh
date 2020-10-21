#!/bin/bash
# compile_run.sh [ TIMES ]

if [ $# == 0 ]; then
  TIMES=1
else
  TIMES="$1"
fi

g++ -fopenmp main.cpp -o main.out

for i in $(seq $TIMES)
do 
  exec ./main.out | tee -a out.log
  echo "-------------------" | tee -a out.log
done

echo "===================" >> out.log