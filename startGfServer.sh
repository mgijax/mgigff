#!/usr/bin/sh

#
#This script starts gfServer 
# Only run the script if the gfServer is not running
#
# To check: ps -ef| grep gfServer
#
cd `dirname $0`
LOG=`pwd`/startGfServer.log
rm -f ${LOG}
touch ${LOG}

working_dir=`pwd`
if [ ! -f Configuration ]
then
  echo "The main configuration file $working_dir/Configuration does not exist"
  exit 1
fi 
. Configuration

date | tee -a ${LOG}


#SERVER_NAME=`uname -n`
SERVER_NAME=$BLAT_HOST
GF_SERVER=`which gfServer`
PORT=$BLAT_PORT

NIBDIR=/data/research/dna/mouse_build_38_nib
LOG_FILE=/data/loads/mgi/mgigff/logs/gfServer.log
rm -f ${LOG_FILE}
touch ${LOG_FILE}

# the set of chromosome nib files to load
CHR="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y MT"

if echo "$GF_SERVER" | grep "not found">/dev/null; then
  echo "gfServer Not Installed on $SERVER_NAME" | tee -a ${LOG}
  exit 1
fi

if [ ! -d ${NIBDIR} ]
then
  echo "${NIBDIR} is not a directory" | tee -a ${LOG}
  exit 1
fi

INPUTFILES=""
# create the string of input files
for c in ${CHR}
do
   INPUTFILES="${INPUTFILES} ${NIBDIR}/chr${c}.nib"
done

echo "gfServer is installed on $SERVER_NAME : $GF_SERVER" | tee -a ${LOG}
echo "Starting blat server: $GF_SERVER start $SERVER_NAME $PORT ${INPUTFILES} -log=$LOG_FILE &" | tee -a ${LOG}
$GF_SERVER start $SERVER_NAME $PORT ${INPUTFILES} -log=$LOG_FILE &

if [ $? -ne 0 ]
then
  echo "gfServer $GF_SERVER failed to start on $SERVER_NAME" | tee -a ${LOG}
  exit 1
fi

echo "gfServer $GF_SERVER started on $SERVER_NAME" | tee -a ${LOG}
exit 0
