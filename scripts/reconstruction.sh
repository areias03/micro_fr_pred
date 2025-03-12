#!/bin/bash -l
#PBS -N mag_reconstruction
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=32:mem=250g
#PBS -q microbiome

eval "$(conda shell.bash hook)"
conda activate carveme
module load cplex

assembly_path="/work/microbiome/global_data_spire/SPIRE/studies/Lloyd-Price_2019_HMP2IBD/psa_megahit/assemblies"
mkdir -p "${assembly_path}/eggnog_mapper"
mkdir -p "${assembly_path}/reconstructions"

for assembly in $assembly_path/*.fa.gz; do
  filename=$(echo "${assembly##*/}")
  echo "${filename}" >> "${assembly_path}/list_of_assemblies.txt"
  IFS='-' read -ra split <<< "$filename"
  sample=$(echo "${split[0]}")
  echo "${sample}" >> "${assembly_path}/list_of_samples.txt"
  sample_list=$(echo "${assembly_path}/list_of_samples.txt")
  # Download EggNOG-mapper data from SPIRE and move to sub-folder
  eggnog_data=$(wget -O "${assembly_path}/eggnog_mapper/${sample}.spire.emapper_annotations.gz" "https://spire.embl.de/download_eggnog/${sample}") 
  tar -xvzf $eggnog_data 
done

# Reconstruction process
cat $sample_list | parallel --dry-run -j 0 carveme --dna "${assembly_path}/{}-assembled.fa.gz" --output "${assembly_path}/reconstructions/{}.xml" --egg "${assembly_path}/eggnog_mapper/{}.spire.emapper_annotations" 
