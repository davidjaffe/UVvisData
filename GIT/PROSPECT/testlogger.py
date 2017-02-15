#!/usr/bin/env python
import Logger,time,sys
if 1: 
    lfn = 'testlog.log'
    sys.stdout = Logger.Logger(fn=lfn)
    print 'output to',lfn
for i in range(10):
    print '\r',i,
    sys.stdout.flush()
    time.sleep(1)
print 'DONE'
