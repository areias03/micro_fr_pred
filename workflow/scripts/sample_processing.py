import os
import os.path as path
import subprocess
import urllib

from jug import TaskGenerator, bvalue

from spirepy import Study

study_name = "Lloyd-Price_2019_HMP2IBD"
data_folder = "/work/microbiome/users/areiasca/micro_fr_pred/data/"
# data_folder = "/home/alexandre/Documents/projects/micro_fr_pred/data/"


# @TaskGenerator
# def process_sample(sample: Sample):
#     print("Sample:\t", sample.id)
#     manifest = sample.manifest
#     manifest.to_csv(
#         path.join(data_folder, study_name, sample.id, "sample_manifest.csv")
#     )
#     return manifest


@TaskGenerator
def download_mag(mag, mag_folder):
    os.makedirs(mag_folder, exist_ok=True)
    print("Started downloading MAG:\t", mag)
    mag_out = path.join(mag_folder, f"{mag}.fa.gz")
    urllib.request.urlretrieve(
        f"https://spire.embl.de/download_file/{mag}",
        mag_out,
    )
    return mag_out


@TaskGenerator
def run_carveme(mag, input, reconstruction_folder):
    print("Started reconstructing MAG:\t", mag)
    os.makedirs(reconstruction_folder, exist_ok=True)
    output = path.join(reconstruction_folder, f"{mag}.xml")
    command = f"carve --dna {input} --output {output} -i M9 -g M9 -v"
    subprocess.check_call(command, shell=True)
    print("Finished processing MAG:\t", mag)


# @TaskGenerator
# def generate_study_manifest(study: Study):
#     list_manifests = []
#     for s in study.samples:
#         list_manifests.append(s.manifest)
#     manifest = pd.concat(list_manifests)
#     manifest.to_csv(path.join(study.folder, "study_manifest.csv"))
#     return manifest


@TaskGenerator
def load_study(study_name: str):
    study = Study(study_name)
    os.makedirs(path.join(data_folder, study_name), exist_ok=True)
    print(f"Loaded study: {study.name}")
    return study


study = load_study(study_name)

for s in bvalue(study).samples:
    print(f"Started Sample: {s}")
    mag_folder = path.join(data_folder, study_name, s.id, "mags")
    reconstruction_folder = path.join(data_folder, study_name, s.id, "reconstructions")
    for mag in bvalue(s.mags["spire_id"].to_list()):
        magf = download_mag(mag, mag_folder)
        run_carveme(mag, magf, reconstruction_folder)

# study_manifest = generate_study_manifest(study)
