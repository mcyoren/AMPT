#!/bin/bash

PW=$(pwd)

cp -r ../AMPT ../AMPT$1 & export PID_C=$!
wait $PID_C
cd ../AMPT$1 || exit
make
sleep 5
./exec & export PID_N=$!
wait $PID_N
cp ana/ampt.dat $PW/output_taxi/ampt_$1.dat & export PID_W=$!
wait $PID_W
cd ../AMPT
rm -r ../AMPT$1
