import os
import os.path as path
import subprocess

from jug import TaskGenerator, bvalue

from spirepy import Study

study_name = "Lloyd-Price_2019_HMP2IBD"
data_folder = "/work/microbiome/users/areiasca/micro_fr_pred/data/"


@TaskGenerator
def load_study(name: str):
    return Study(name)


@TaskGenerator
def reconstruct_mag(mag: str, mag_folder: str, reconstruction_folder: str):
    os.makedirs(reconstruction_folder, exist_ok=True)
    print(f"Started reconstructing:\t {mag}")
    input = path.join(mag_folder, f"{mag}.fa.gz")
    output = path.join(reconstruction_folder, f"{mag}.xml")
    command = f"carve --dna {input} --output {output} -i M9 -g M9 -v"
    subprocess.check_call(command, shell=True)
    print(f"Finished reconstructing:\t {mag}")


study = load_study(study_name)

for m in bvalue(study).mags.iter_rows(named=True):
    reconstruct_mag(
        m["genome_id"],
        path.join(data_folder, study_name, m["derived_from_sample"], "mags"),
        path.join(data_folder, study_name, m["derived_from_sample"], "reconstructions"),
    )
