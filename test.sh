#!/bin/bash

python construct_protein_graph.py test.pdb > test.graph

function check {
  if [ $1 != $3 ]; then
      echo "Test Failed"
      echo "$1 $2 $3 $4"
      exit
  fi
  if [ $2 != $4 ]; then
      echo "Test Failed"
      echo "$1 $2 $3 $4"
      exit
  fi
}

CORRECT=("2 3" "2 4" "2 5" "3 2" "4 2" "5 2" "6 7" "6 8" "7 6" "7 8" "8 6" "8 7")
for i in `seq 12`
  do
  OUT=`head -n $i test.graph | tail -n 1`
  IDX=`expr $i - 1`
  check $OUT ${CORRECT[$IDX]}
done
echo "Test Success"
