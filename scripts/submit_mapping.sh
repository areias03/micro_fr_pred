#!/bin/bash -l
#PBS -N mag_spire_mapping
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=16:mem=125g
#PBS -q microbiome

study_name="Valles-Colomer_2023_transmission"
project_path="/work/microbiome/users/areiasca/functional_prediction"

mkdir -p "$project_path/data/$study_name/"
wget "https://spire.embl.de/api/study/$study_name?format=csv" -O "$project_path/data/$study_name.csv"

output="$project_path/data/$study_name/"
study_data="$project_path/data/$study_name.csv"

python "$project_path/scripts/study_mag_mapping.py" "$study_data"

list_of_mags="$project_path/data/list_of_targets.txt"

while IFS=",", read -r mag sample
do
    target=$(echo "${output}${sample}");
    mkdir -p $target;
    cp $mag $target &
done < "$list_of_mags"

# cat $list_of_mags | parallel --dry-run -j 300% cp {} $output
