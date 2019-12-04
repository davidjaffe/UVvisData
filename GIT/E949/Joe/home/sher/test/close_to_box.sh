#!/bin/tcsh
#
#    File:  close_to_box.sh
# Purpose:  find events close to the box
#   Usage:  > close_to_box.sh {background_type} {year}
#              e.g.  close_to_box.sh kp2 98
#
# History:  2001 Aug  PB  created
#

#
#  Command line syntax.
#
if ($#argv < 2) then
   echo "\n Usage:  close_to_box.sh {background_type} {year} \n"
   exit 1
endif
set type = $1
set year = $2
if ($type == "kp2") then
   set str1 = skim1
   set str2 = skim4
else if ($type == "km2") then
   set str1 = skim2
   set str2 = skim5
else if ($type == "bm") then
   set str1 = skim3
   set str2 = skim6
endif

#
#  Run paw job.
#
set logfile = close_to_box.$type.$year.log
echo "\n `date` \n" > $logfile

set machine = `hostname`
if (`echo $machine | grep -i sitka | wc -w`) then
   source $E787_SETUP_FILE
   setup pass2 98
endif

set command = `which paw | awk '{print $NF}'`
nice +2 $command <<+ >> $logfile
0
0
exec chain $str1 $year
exec chain $str2 $year
chain data $str1 $str2
cd //data
nt/loop 1 close_to_box.f(0.)
exit
+

\mv close_to_box.dat close_to_box.$type.$year.dat

echo "\n `date` \n" >> $logfile

exit 0

