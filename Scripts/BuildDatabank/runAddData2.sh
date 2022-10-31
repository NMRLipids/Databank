#!/bin/bash
#Loop AddData.py over many simulations 

for i in {677..741}
do
    echo $i
    python AddData.py -f ./info_files/$i/info.yaml
done
