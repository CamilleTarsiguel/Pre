#!/bin/sh

grepstr='12';
if [[ $# -eq 1 ]]; then
	grepstr=$1;
fi
for dir in ./results/$grepstr*/
do
  dir=${dir%/}
  dir=${dir##*/}
#   echo $dir
#  cat results/$dir/bestres.txt
  ll=`tail results/$dir/bestres.txt -n 1`
  mota=${ll:54:4}
  motp=${ll:60:4}
  mt=${ll:23:3}
  ml=${ll:30:3}
  id=${ll:44:4}
  fm=${ll:48:3}
  rcl=${ll:1:4}
  prc=${ll:7:4}

  echo $dir  $mota  $motp  $mt  $ml  $id  $fm  $rcl  $prc
done
