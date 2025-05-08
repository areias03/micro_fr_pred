import os
import os.path as path
import subprocess
import urllib

import pandas as pd
from jug import TaskGenerator, bvalue

from micro_fr_pred import Sample, Study

study_name = "Lloyd-Price_2019_HMP2IBD"
data_folder = "/work/microbiome/users/areiasca/micro_fr_pred/data/"
# data_folder = "/home/alexandre/Documents/projects/micro_fr_pred/data/"


@TaskGenerator
def process_sample(sample: Sample):
    print("Sample:\t", sample.id)
    manifest = sample.manifest
    manifest.to_csv(path.join(s.out_folder, "sample_manifest.csv"))
    return manifest


@TaskGenerator
def download_mag(mag, mag_folder):
    os.makedirs(mag_folder, exist_ok=True)
    print("Started MAG:\t", mag)
    mag_out = path.join(mag_folder, f"{mag}.fa.gz")
    urllib.request.urlretrieve(
        f"https://spire.embl.de/download_file/{mag}",
        mag_out,
    )
    return mag_out


@TaskGenerator
def run_carveme(mag, input, reconstruction_folder):
    os.makedirs(reconstruction_folder, exist_ok=True)
    output = path.join(reconstruction_folder, f"{mag}.xml")
    command = f"carve --dna {input} --output {output} -i M9 -g M9 -v"
    subprocess.check_call(command, shell=True)
    print("Finished processing MAG:\t", mag)


@TaskGenerator
def generate_study_manifest(study: Study):
    list_manifests = []
    for s in study.samples:
        list_manifests.append(s.manifest)
    manifest = pd.concat(list_manifests)
    manifest.to_csv(path.join(study.folder, "study_manifest.csv"))
    return manifest


@TaskGenerator
def load_study(study_name: str):
    study = Study(study_name, path.join(data_folder, study_name))
    return study


study = load_study(study_name)

for s in bvalue(study).samples:
    manifest = process_sample(s)
    mag_folder = path.join(s.out_folder, "mags")
    reconstruction_folder = path.join(s.out_folder, "reconstructions")
    for mag in bvalue(s.mags["genome_id"].tolist()):
        magf = download_mag(mag, mag_folder)
        run_carveme(mag, magf, reconstruction_folder)
    print(manifest)

study_manifest = generate_study_manifest(study)
