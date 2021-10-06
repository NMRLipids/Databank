#!/bin/bash
#Loop AddData.py over many simulations 

for i in {179..183}
do
    python3.6 AddData.py -f ./info_files/$i/*.yaml
done
