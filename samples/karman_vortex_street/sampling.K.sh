#!/bin/sh
# [ shell script redraw-step.sh; update: 2022-03-31 ]
# call this script from directory ${WRKDIR}
# 
#BLDDIR=${HOME}/danse/danse-1.0r3/danse-1.0r3
#BLDDIR=../
BLDDIR=../current/danse-1.0r3/
BINDIR=${HOME}/bin
#WRKDIR=${HOME}/tmp/step-data
WRKDIR=./
MODDIR=../current/danse-1.0r3/models
MODEL="mod_step"
POSTERINPUT="posterinput.K"
GNUINPUT="gnuinput"

${BINDIR}/solver.do -b

#NUMBER="${1}"
NUMBER="-1"
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
                sleep 5
	done
}

if [ -x gnuplot ] \
|| [ -x /bin/gnuplot ] \
|| [ -x /usr/bin/gnuplot ] \
|| [ -x /usr/local/bin/gnuplot ]
then
	if [ -x ${BINDIR}/poster.do ]
	then
		iterate
	else
		echo "program poster not available !!!"
		echo "trying to compile ... "
		cd ${BLDDIR}
		rm model.c
		ln -s ${MODDIR}${MODEL} model.c
        	make clean
        	make poster

		if [ -x ${BINDIR}/poster.do ]
		then
			cd ${WRKDIR}
			iterate
		else
			echo "program poster not available !!!"
			exit 0
		fi
	fi
else
	echo "program gnuplot not available !!!"
	exit 0
fi
