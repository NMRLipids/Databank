#!/bin/bash
# Check if all systems in the NMRlipids databank contain an unique ID
# 
# If the argument "generate" is provided, README.yaml files where
# indices are duplicated or missing will be modified to correct it.
#
# Contact:
#  Fabián Suárez-Lestón
#  fabian.suarez.leston@usc.es

cd Data/Simulations

paths=()
ids=()

# Iterate over all systems
for file in `find . -name "README.yaml"`; do
 # The path to the file
 path=$( echo ${file} | rev | cut -d"/" -f2- | rev )
 paths+=(${path})
 # The ID of the file
 file_ID=$( grep '^ID:' ${file} | cut -d" " -f2 )
 # If the file does not contain an ID
 if ! [ $file_ID ]; then
  ids+=(0)
 else
  ids+=(${file_ID})
 fi
done

# The list on unique IDs in the databank
idlist=($(printf '%s\n' "${ids[@]}" | sort -u | sort -n ))

# The last ID used
max_id=$( cat COUNTER_ID )

# Find duplicates in the array of IDs; from https://stackoverflow.com/questions/22055238/search-duplicate-element-array
duplicates=($( printf '%s\n' "${ids[@]}"|awk '!($0 in seen){seen[$0];next} 1' ))

if [[ ${#duplicates[@]} -eq 0 ]] && [[ ${idlist[0]} -ne 0 ]]; then
 echo "### :white_check_mark: All systems have an unique ID associated" >> $GITHUB_STEP_SUMMARY

 # No IDs were generated
 echo "newids=false" >> $GITHUB_OUTPUT

else
 echo "### :warning: Duplicates and/or missing IDs have been found" >> $GITHUB_STEP_SUMMARY
 echo " " >> $GITHUB_STEP_SUMMARY

 fixlist=()
 unique_ids=()

 # Repeated IDs
 if [[ ${#duplicates[@]} -ne 0 ]]; then
  unique_ids=($( printf "%s\n" "${duplicates[@]}" | sort -u ))
 fi

 # No ID assigned
 if [[ ${idlist[0]} -eq 0 ]] && ( [[ ${unique_ids[0]} -ne 0 ]] || ! [[ ${unique_ids[0]} ]] ) ; then
  unique_ids=(0 ${unique_ids[@]})
 fi

 # Find all duplicated and missing IDs; ID matches the index of the array
 for j in $( seq ${#paths[@]} ); do
  for u in ${unique_ids[@]}; do
   if [[ ${ids[$((${j}-1))]} -eq ${u} ]]; then
    fixlist[${u}]+=${paths[$((${j}-1))]}" "
    break
   fi
  done
 done

 # Print the results
 counter=0
 for u in $( printf '%s\n' "${unique_ids[@]}" | sort -r  ); do
  list=(${fixlist[${u}]})
  if [[ ${u} -eq 0 ]]; then
   echo "#### Systems with no ID assigned:" >> $GITHUB_STEP_SUMMARY
   for item in ${list[@]}; do
    echo \`${item}\` >> $GITHUB_STEP_SUMMARY
   done
   counter=$((${counter}+${#list[@]}))
  else
   echo "#### Systems with duplicated index ${u}": 
   for item in ${list[@]}; do
    echo \`${item}\` >> $GITHUB_STEP_SUMMARY
   done
   counter=$((${counter}+${#list[@]}-1))
  fi
 done
 echo "#### Number of systems to be fixed: "${counter} >> $GITHUB_STEP_SUMMARY
 echo " " >> $GITHUB_STEP_SUMMARY

 if [[ "${1}" = "generate" ]]; then
  # Generate a list of valid IDs
  new_id=( $( seq $((${max_id}+1)) $((${max_id}+${counter})) ) )
  counter=0
  # Apply changes in the IDs of the systems
  for u in $( printf '%s\n' "${unique_ids[@]}" | sort -r  ); do
   list=(${fixlist[${u}]})
   if [[ ${u} -eq 0 ]]; then
    echo "#### Fixing missing indices..." >> $GITHUB_STEP_SUMMARY
    ordered=($( for path in ${list[@]}; do echo ${path} $( date -d $( grep "DATEOFRUNNING:" ${path}/README.yaml | cut -d" " -f2 | awk -F' |/' '{printf "%s-%s-%s %s",$3,$2,$1,$4}' ) +%s ) $( git log -1 --pretty="format:%ct" "${path}/README.yaml" ); done | sort -n -k2,3 | cut -d" " -f1 ))
    for path in ${ordered[@]}; do
     echo "#### ID "${new_id[counter]}" :arrow_right: \`"${path}"\`" >> $GITHUB_STEP_SUMMARY
     echo "ID: "${new_id[counter]} >> ${path}/README.yaml
     counter=$((${counter}+1))
    done
   else
    echo "#### Fixing duplicated index ${u}..." >> $GITHUB_STEP_SUMMARY
    ordered=($( for path in ${list[@]}; do echo ${path} $( date -d $( grep "DATEOFRUNNING:" ${path}/README.yaml | cut -d" " -f2 | awk -F' |/' '{printf "%s-%s-%s %s",$3,$2,$1,$4}' ) +%s ) $( git log -1 --pretty="format:%ct" "${path}/README.yaml" ); done | sort -n -k2,3 | cut -d" " -f1 ))
    for k in $(seq 1 $((${#ordered[@]}-1))); do
     echo "#### ID "${new_id[counter]}" :arrow_right: \`"${ordered[k]}"\`" >> $GITHUB_STEP_SUMMARY
     sed -i"" "s/ID:.*/ID: ${new_id[counter]}/" ${ordered[k]}/README.yaml
     counter=$((${counter}+1))
    done
   fi
  done
 
  # Update the list of used IDs
  echo ${new_id[$((${counter}-1))]} > COUNTER_ID

  # New IDs were generated
  echo "newids=true" >> $GITHUB_OUTPUT
 fi
fi
