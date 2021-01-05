#!/bin/bash
for i in {0..19}
do
  echo $i
  ncl data-extract_1_tsmpBonnRadar.ncl expt=$1 inst=$i
done
exit
