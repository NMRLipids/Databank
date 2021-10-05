#!/bin/bash
#Loop AddData.py over many simulations 

for f in info_files/3/*.yaml
do
    /home/akiirikk/anaconda3/bin/python3.7 AddData.py -f $f
done
