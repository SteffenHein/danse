#!/bin/sh
# [ shell script redraw-step.sh; update: 2022-03-31 ]
# call this script from directory ${WRKDIR}
# 
#BLDDIR=${HOME}/tmp/danse/danse-1.0r3/danse-1.0r3
BLDDIR=../
BINDIR=${HOME}/bin
#WRKDIR=${HOME}/tmp/step-data
WRKDIR=../work
MODDIR=../models
MODEL="mod_step"
POSTERINPUT="posterinput.B"
GNUINPUT="gnuinput"

${BINDIR}/solver.do -b

#NUMBER="${1}"
NUMBER="0"
iterate( )
{
	while [ "${NUMBER}" != 1000 ]; do
		NUMBER=`expr ${NUMBER} + 1`

		${BINDIR}/poster.do < ${POSTERINPUT}
#cp flw.dsp flw.dsp${NUMBER}
		cp flw.dat ${NUMBER}"_TU"
		gnuplot < ${GNUINPUT}
 		cp flw.pdf ${NUMBER}"_TU.pdf"

		if [ ! -f ${WRKDIR}/dsc.pid ]
		then
			echo "Terminated"
			exit 0
		fi
                sleep .2
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
