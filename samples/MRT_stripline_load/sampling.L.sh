#!/bin/sh
# [ shell script redraw-step.sh; update: 2022-03-31 ]
# call this script from directory ${WRKDIR}
# 
#BLDDIR=${HOME}/danse/danse-1.0r3/danse-1.0r3
BLDDIR=../
BINDIR=${HOME}/bin
#WRKDIR=${HOME}/tmp/load-data
WRKDIR=./
MODDIR=../
MODEL="model.c"
POSTERINPUT="posterinput.L"
GNUINPUT="gnuinput"

${BINDIR}/solver.do -b

#NUMBER="${1}"
NUMBER="0"
iterate( )
{
	while [ "${NUMBER}" != 10000 ]; do
		NUMBER=`expr ${NUMBER} + 1`

		${BINDIR}/poster.do < ${POSTERINPUT}
#cp flw.dsp flw.dsp${NUMBER}
		cp flw.dat ${NUMBER}"_TU"
		gnuplot < ${GNUINPUT}
 		ps2pdf flw.ps ${NUMBER}"_TU.pdf"

		if [ ! -f dsc.pid ]
		then
			echo "Terminated"
			exit 0
		fi
                sleep 0.1
	done
}
iterate
exit 0
