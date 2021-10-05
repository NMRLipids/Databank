#!/bin/bash
#Loop AddData.py over many simulations 

for i in {7..137}
do
    python3.6 AddData.py -f ./info_files/$i/*.yaml
done
