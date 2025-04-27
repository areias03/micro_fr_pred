import os.path as path
import pandas as pd
from jug import TaskGenerator
from micro_fr_pred import Study, Sample

study_name = "Lloyd-Price_2019_HMP2IBD"
# data_folder = "/work/microbiome/users/areiasca/micro_fr_pred/data/"
data_folder = "~/microbiome/users/areiasca/micro_fr_pred/data/"


@TaskGenerator
def process_sample(sample: Sample):
    sample.download_mags()
    sample.reconstructions
    manifest = sample.manifest
    return manifest


@TaskGenerator
def generate_study_manifest(study: Study):
    list_manifests = []
    for s in study.samples:
        list_manifests.append(s)
    manifest = pd.concat(list_manifests)
    return manifest


@TaskGenerator
def main():
    study = Study(study_name, path.join(data_folder, study_name))
    for s in study.samples:
        sample_manifest = process_sample(s)
        sample_manifest.to_csv(path.join(s.out_folder, "sample_manifest.csv"))
    study_manifest = generate_study_manifest(study)
    study_manifest.to_csv(path.join(study.folder, "study_manifest.csv"))


if __name__ == "__main__":
    main()
