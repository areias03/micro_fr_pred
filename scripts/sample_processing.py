import os.path as path
import pandas as pd
from jug import TaskGenerator
from jug import bvalue
from micro_fr_pred import Study, Sample

study_name = "Lloyd-Price_2019_HMP2IBD"
# data_folder = "/work/microbiome/users/areiasca/micro_fr_pred/data/"
data_folder = "/home/alexandre/Documents/projects/micro_fr_pred/data/"


@TaskGenerator
def process_sample(sample: Sample):
    print(sample.id)
    sample.download_mags()
    sample.reconstructions
    manifest = sample.manifest
    manifest.to_csv(path.join(s.out_folder, "sample_manifest.csv"))
    return manifest


@TaskGenerator
def generate_study_manifest(study: Study):
    list_manifests = []
    for s in study.samples:
        list_manifests.append(s)
    manifest = pd.concat(list_manifests)
    manifest.to_csv(path.join(study.folder, "study_manifest.csv"))
    return manifest


@TaskGenerator
def load_study(study_name: str):
    study = Study(study_name, path.join(data_folder, study_name))
    return study


study = load_study(study_name)

for s in bvalue(study).samples:
    sample_manifest = process_sample(s)
study_manifest = generate_study_manifest(study)
