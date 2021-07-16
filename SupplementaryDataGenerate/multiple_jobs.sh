#! /usr/bin/env bash

set -o monitor 
# means: run background processes in a separate processes...

trap add_next_job CHLD 
# execute add_next_job when we receive a child complete signal

todo_array=($(find "/mnt/ResearchData/01_SCOPUS/00_SCOPUS_20190421/ANI-ITEM/" -name "*.zip")) # places output into an array

unzip -p /mnt/ResearchData/01_SCOPUS/00_SCOPUS_20190421/ANI-ITEM/2019-4-21_FULLFORMAT_010-xml-1.zip | 

index=0
max_jobs=64

function add_next_job {
    # if still jobs to do then add one
    if [[ $index -lt ${#todo_array[*]} ]]
    # apparently stackoverflow doesn't like bash syntax
    # the hash in the if is not a comment - rather it's bash awkward way of getting its length
    then
        echo adding job ${todo_array[$index]}
        do_job ${todo_array[$index]} & 
        # replace the line above with the command you want
        index=$(($index+1))
    fi
}

function do_job {
#    echo "starting job $1.out"
    echo ${1##*/}
    unzip -p $1 | pypy scopus_xml_parsing_20191021.py > ./pased_scopus/${1##*/}.out 
    sleep 0.5
}

# add initial set of jobs
while [[ $index -lt $max_jobs ]]
do
    add_next_job
done

# wait for all jobs to complete
wait
echo "done"
