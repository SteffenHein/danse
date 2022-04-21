#!/bin/sh
# on time interval $SLEEP [seconds], displays $LOGDIR/dsc.log
# on console $CONSOLE
SLEEP=60
#LOGDIR=.
LOGDIR=/var/log		#usually the directory, where DSC job is running 
CONSOLE=/dev/ttyv1	#/dev/ttyv0...7
while [ 1 ]; do
	cat $LOGDIR/dsc.log > $CONSOLE
	sleep $SLEEP
done

