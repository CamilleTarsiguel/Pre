#!/bin/sh
cat results/$1/bestres.txt
ll=`tail results/$1/bestres.txt -n 1`
mota=${ll:54:4}
motp=${ll:60:4}
mt=${ll:23:3}
ml=${ll:30:3}
id=${ll:44:4}
fm=${ll:48:3}
rcl=${ll:1:4}
prc=${ll:7:4}

echo $mota \& $motp \& $mt \& $ml \& $id \& $fm \& $rcl \& $prc
