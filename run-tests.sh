#!/bin/bash
cd `dirname $0`
SCRIPT_DIR=`pwd`

if [ ! -e $SCRIPT_DIR/databases ];
then
	./staramr.py db build
fi

python -m unittest discover
