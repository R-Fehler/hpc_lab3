#!/usr/bin/env bash
if [[ "$OSTYPE" == "linux-gnu" ]]; then
        alias gtime=time
fi
mkdir -p log
LOOP=3
NTIMESTEPS=100
num_threads_in_x=(1 2 4) #bash Array
num_threads_in_y=(1 2 2)
ARRAYSIZE_PER_THREAD_X=400
ARRAYSIZE_PER_THREAD_Y=800
PATH_TO="mpiexec -np 4 ./build/gameoflife_mpi"
    echo "FeldNr. Zeit in sek CPU %">log/logtime 
    echo "LOG:" > log/logprogram

    for (( i = 0; i < $LOOP; i++ )); do
    for (( j = 0; j < 3; j++ )); do
    cnt=$(($j+1))
   (gtime -f "$cnt %e %P" $PATH_TO ${num_threads_in_x[$j]} ${num_threads_in_y[$j]} $ARRAYSIZE_PER_THREAD_X $ARRAYSIZE_PER_THREAD_Y $NTIMESTEPS) >>log/logprogram 2>>log/logtime
    done
    done
awk '$1=="1" {printf"%s, %s\n", $2, $3}' log/logtime > log/logtime1.csv
awk '$1=="2" {printf"%s, %s\n", $2, $3}' log/logtime > log/logtime2.csv
awk '$1=="3" {printf"%s, %s\n", $2, $3}' log/logtime > log/logtime3.csv


    

