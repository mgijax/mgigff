#!/bin/sh
cd `dirname $0`
MGIGFFHOME=`pwd`
export MGIGFFHOME
. ${MGIGFFHOME}/Configuration
${PYTHON} ${MGIGFFHOME}/bin/mgigff.py $*
