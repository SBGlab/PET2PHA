#!/bin/bash
export PATH=$PATH:/home/alvaro/Matlab_2022/bin/

nohup matlab -nodisplay -nodesktop -r "run('/home/alvaro/github/PET2PHA/code/utils/run_gcFront.m')" > /home/alvaro/github/PET2PHA/data/gcfront.txt &
