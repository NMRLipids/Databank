#!/bin/bash
# Check if there has been conducted a search for compatible experiments for all
# simulations in the NMRlipids databank
#
# Contact:
#  Fabián Suárez-Lestón
#  fabian.suarez.leston@usc.es

cd Data/Simulations

paths=()

# Iterate over all systems
for file in `find . -name "README.yaml"`; do
 if ! [ $( grep "EXPERIMENT" ${file} ) ]; then
  # The path to the file
  path=$( echo ${file} | rev | cut -d"/" -f2- | rev )
  paths+=(${path})
 fi
done

if [[ ${#paths[@]} -ne 0 ]]; then
 echo "### :warning: There has been no search for compatible experiments for:" >> $GITHUB_STEP_SUMMARY
 for system in ${paths[@]}; do
  echo \`${system}\` >> $GITHUB_STEP_SUMMARY
 done
else
 echo "### :white_check_mark: A search for compatible experiments has been performed for all systems." >> $GITHUB_STEP_SUMMARY
fi
