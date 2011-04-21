#!/bin/sh
echo 'Running minimal solver test'
rm -fr iter*
mpirun -np 1 ./minimal -p parameters -m usquare.xml
rm -fr iter*
mpirun -np 4 ./minimal -p parameters -m usquare.xml
rm -fr iter*
mpirun -np 9 ./minimal -p parameters -m usquare.xml