#!/bin/bash

current=0

if [ -f ".last_run_format" ]; then
  current=`cat .last_run_format`
fi

for file in src/*.cpp src/{algorithm,classes,graph,python,utils}/*.cpp src/algorithm/{electron_assignment,graph}/*.cpp src/python/{classes,graph}/*.cpp include/indigox/{algorithm,classes,graph,python,utils}/*.hpp include/indigox/algorithm/{electron_assignment,graph}/*.hpp
do
last_modified=`stat -c "%Y" $file`
if [ $(($current-$last_modified)) -lt 0 ]; then
    clang-format -i $file
fi
done

date +%s > .last_run_format
