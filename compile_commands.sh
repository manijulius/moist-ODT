!#/bin/bash

#compile and run init program
gfortran -O3 Labinit_moist.f  -o Labinit_moist.exe
./Labinit_moist.exe

#compile the exp program
gfortran -O3 LabExp_moist.f  -o LabExp_moist.exe

#run experiment to steady state
nohup ./LabExp_moist.exe default steady_state 001 &> steady_state.log &
mkdir -p input/steady_state
cp  output/steady_state/001/*   input/steady_state/

#test experiment starting from steady_state output files
nohup ./LabExp_moist.exe steady_state test 001 &> 001.log &
