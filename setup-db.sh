#!/bin/bash
cd `dirname $0`
SCRIPT_DIR=`pwd`
export PYTHONPATH=$SCRIPT_DIR:$PYTHONPATH

if [ ! -e $SCRIPT_DIR/staramr/databases/data/dist ];
then
	python bin/staramr db build
fi
