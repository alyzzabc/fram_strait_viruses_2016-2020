#!/bin/bash

#####
### Run orthofinder 
#####

## Read input

############################################################
# Help                                                     #
############################################################
Help()
{
   # Display Help
   echo "Calculate coverage of virus per Gbp of metagenome reads"
   echo
   echo "Syntax: calc_cvpg.sh [-i|-n|-t]"
   echo "options:"
   echo "-i / --indir     Input directory"
   echo "-t / --threads   Number of threds"
   echo "h     Print this Help"
   echo
}

############################################################
# Main                                                     #
############################################################

SHORT=i:,t:,h
LONG=indir:,threads:,outdir:,help
OPTS=$(getopt --options $SHORT --longoptions $LONG -- "$@")

VALID_ARGUMENTS=$# # Returns the count of arguments that are in short or long options

if [ "$VALID_ARGUMENTS" -eq 0 ]; then
  Help
fi

eval set -- "$OPTS"

while :
do
    case $1 in
        -i | --indir ) 
        indir="$2" 
        shift 2 
        ;;
        -t | --threads ) 
        threads="$2" 
        shift 2 
        ;;
        -h | --help ) 
        Help
        exit 2
        ;;
        --)
        shift;
        break
        ;;
        *)
        echo "Unexpected option: $1"
        Help
        exit 2
        ;;
  esac
done

cd ${indir}

echo "checking that the necessary files exist"
count1=$( ls -l *.bam | wc -l )

if [[ $count1 != 0 && -s FRAM_RAS_name_matching_and_size.txt && -s calc_coverage_statistics.py ]];
then
    echo "Files found. Continuing with coverage statistic estimation"
    for i in *.bam ; do j=$( echo $i | sed 's/.bam//' ) ; samtools depth -a -@ ${threads} "$i" > "$j".bed ; done
else
    echo "No bam files were found in input directory"
fi

count2=$( ls -l *.bed | wc -l )
if [[ $count2 != 0 ]];
then
    for i in *.bed; do j=$( echo $i | sed 's/.bed//' ) ; python calc_coverage_statistics.py -i "$i" -o "$j"_coverage_statistics.txt ; done
else
    echo "Bam to bedgraph conversion failed"
fi

count3=$( ls -l *.bed | wc -l )
if [[ $count3 != 0 ]];
then
    grep -v 'Sequence' *coverage_statistics.txt | sed 's/.filter.*txt:/\t/' | \
    awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$5]=$1"\t"$3"\t"$4;next}{if ($1 in a) print $1,$2,$3,$4,a[$1]}' FRAM_RAS_name_matching_and_size.txt - | \
    awk 'BEGIN{FS=OFS="\t"}{print $1,$6,$2,$3,$4,$5,$5/$7,$5/$8}' | sed '1i Library_name\tRAS_id\tRead_name\tLength\tBreadth_coverage\tMean_depth_coverage\tCVPG\tCV_per_genome' > FRAM_RAS_read_CVPG.txt
else
    echo "Coverage statistics from bedgraph files was not completed successfully"
fi
