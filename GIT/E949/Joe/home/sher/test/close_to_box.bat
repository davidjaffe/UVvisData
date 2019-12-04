#!/bin/tcsh

foreach sample (98)
   echo "submitting close_to_box km2 $sample ..."
   close_to_box.sh km2 $sample >& /dev/null
   echo "submitting close_to_box kp2 $sample ..."
   close_to_box.sh kp2 $sample >& /dev/null
   echo "submitting close_to_box bm $sample ..."
   close_to_box.sh bm $sample >& /dev/null
end

exit 0

