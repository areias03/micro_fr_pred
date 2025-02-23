#!/bin/bash -l
#PBS -N mag_reconstruction
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=31:mem=248g
#PBS -q microbiome

assembly_path="/work/microbiome/global_data_spire/SPIRE/studies/Lloyd-Price_2019_HMP2IBD/psa_megahit/assemblies"
mkdir -p "${assembly_path}/eggnog_mapper"
mkdir -p "${assembly_path}/reconstructions"

for assembly in $assembly_path/*.fa.gz; do
  filename=$(echo "${assembly##*/}")
  IFS='-' read -ra split <<< "$filename"
  sample=$(echo "${split[0]}")
  # Download EggNOG-mapper data from SPIRE and move to sub-folder
  eggnog_data=$(wget "https://spire.embl.de/download_eggnog/${sample}" -O "${assembly_path}/eggnog_mapper/${sample}.spire.emapper_annotations") &
  echo "${eggnog_data}"
  # tar -xzf $eggnog_data &
  # Reconstruction process
  conda activate carveme
  module load cplex
  # parallel --dry-run -j 0 carveme $assembly --output "${assembly_path}/reconstructions/${sample}.xml" --egg "${assembly_path}/eggnog_mapper/${sample}.spire.emapper_annotations" 
done

