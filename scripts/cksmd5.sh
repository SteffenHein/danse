#!/bin/sh
CHECKSUM=$1
RMFILES=${CHECKSUM}" .md5"
#
cd ../models
rm -f ${RMFILES}
md5sum * > ${CHECKSUM}
cd ../bin
rm -f ${RMFILES}
md5sum * > ${CHECKSUM}
cd ../doc
rm -f ${RMFILES}
md5sum * > ${CHECKSUM}
cd ../former
rm -f ${RMFILES}
md5sum * > ${CHECKSUM}
cd ../poster
rm -f ${RMFILES}
md5sum * > ${CHECKSUM}
cd ../solver
rm -f ${RMFILES}
md5sum * > ${CHECKSUM}
cd ../expm
rm -f ${RMFILES}
md5sum * > ${CHECKSUM}
cd ../math
rm -f ${RMFILES}
md5sum * > ${CHECKSUM}
cd ../objects
rm -f ${RMFILES}
md5sum * > ${CHECKSUM}
cd ../history
rm -f ${RMFILES}
md5sum * > ${CHECKSUM}
cd ../samples
rm -f ${RMFILES}
md5sum * > ${CHECKSUM}
cd ../scripts
rm -f ${RMFILES}
md5sum * > ${CHECKSUM}
cd ../tools
rm -f ${RMFILES}
md5sum * > ${CHECKSUM}
cd ../
rm -f ${RMFILES}
md5sum * > ${CHECKSUM}
