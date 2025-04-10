import os
from micom.viz import plot_growth
from micom.workflows import grow

from micom.qiime_formats import load_qiime_medium
from micom.workflows import build, build_database

from micro_fr_pred import Study
from micro_fr_pred.spire import genome_metadata

data_folder = "/home/alexandre/Documents/projects/micro_fr_pred/data"
study_name = "Lloyd-Price_2019_HMP2IBD"

study = Study(study_name, os.path.join(data_folder, study_name))
sample = study.samples[-1]

# Tests

print("###Study ###")
print(f"Name: {study.name}\n")
print(f"Metadata:\n {study.metadata}\n")

print("###Sample ###")
print(f"Name: {sample.id}\n")
print(f"MAGs:\n {sample.mags}\n")
print(f"Abundances:\n {sample.get_abundances()}\n")
print(f"Manifest:\n {sample.manifest.file[0]}\n")

sample.manifest.to_csv(f"{sample.out_folder}sample_manifest.csv")
