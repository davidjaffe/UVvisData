#!/bin/tcsh

set id_list = (9597inside 9597insidelo 9597insidehi)
set id_list = ($id_list 98inside 98insidelo 98insidehi 98outside 98outside0 98total)
set id_list = ($id_list 9598inside 9598insidelo 9598insidehi 9598outside 9598outside0 9598total)

foreach id ($id_list)
  
  wait_here:
  @ nrunning = `ps -ef | grep ecl | grep exe | grep -v grep | wc -l`
  if ($nrunning < 4) then
    echo "running $id..."  
    ecl$id.exe >& ecl$id.dat &
    sleep 2
    goto next_job
  endif
  
  sleep 30
  goto wait_here

  next_job:
end
      
exit
