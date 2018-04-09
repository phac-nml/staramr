#!/bin/bash
cd `dirname $0`
SCRIPT_DIR=`pwd`

if [ ! -e $SCRIPT_DIR/databases ];
then
	python bin/staramr db build
fi

nosetests
