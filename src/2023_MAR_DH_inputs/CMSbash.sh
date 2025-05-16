#!/bin/bash
for filename in job*; do
  [ -e "$filename" ] || continue
  qsub "$filename"
  echo "$filename has been submitted"
  i=$(qselect -u holstein | wc -l)
  echo "The number of submitted jobs is: $i"
  while [ $i -ge 4 ];
  do
    echo "Wait..."
    read -t 120
    i=$(qselect -u holstein | wc -l)
    echo "The number of submitted jobs is: $i"
  done
  #echo "submit job $filename"
done
