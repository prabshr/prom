#!/bin/bash
#Main code
rundir=$PROM/scratch
for inst in {0..19}
do
  echo "Processing instance "  $inst
  cp cluma2job $rundir/cluma2job_$inst
  sed "s,__inst__,$inst," -i ${rundir}/cluma2job_$inst
  sed "s,__inst__,$inst," -i ${rundir}/cluma2job_$inst
  qsub $rundir/cluma2job_$inst
done

exit 0
