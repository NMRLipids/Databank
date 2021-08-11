awk '{if(NR != 2 && NR != 3) print $0}' /media/akiirikk/DATADRIVE1/tietokanta/NMRlipidsIVPEandPG/scripts/DataBankINFO/POPScharmmINFO.py > info.yaml
perl -i -pe 's/@//g' info.yaml
perl -i -pe 's/=/: /g' info.yaml
perl -i -pe 's/"//g' info.yaml
perl -i -pe 's/dir_wrk/DIR_WRK/g' info.yaml
perl -i -nle 'print if !/NPOPC/' info.yaml
perl -i -nle 'print if !/NPOPG/' info.yaml
perl -i -nle 'print if !/NPOPS/' info.yaml
perl -i -nle 'print if !/NPOPE/' info.yaml
perl -i -nle 'print if !/NPOT/' info.yaml
perl -i -nle 'print if !/NSOD/' info.yaml
perl -i -nle 'print if !/NCLA/' info.yaml
perl -i -nle 'print if !/NCAL/' info.yaml
perl -i -nle 'print if !/NSOL/' info.yaml
perl -i -nle 'print if !/TEMPERATURE/' info.yaml
perl -i -nle 'print if !/TRJLENGTH/' info.yaml
perl -i -nle 'print if !/#NMRLIPIDS/' info.yaml
perl -i -nle 'print if !/user_information/' info.yaml
perl -i -nle 'print if !/SIM/' info.yaml
perl -i -pe 's/,/: /g' info.yaml
awk '{if($0 ~ "MAPPING") print "MAPPING_DICT:\n  " $2; else print $0 }' info.yaml > tmp.yaml
mv tmp.yaml info.yaml
