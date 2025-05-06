import os
import os.path as path
import subprocess

from jug import TaskGenerator, bvalue

from micro_fr_pred import Study

study_name = "Lloyd-Price_2019_HMP2IBD"
data_folder = "/work/microbiome/users/areiasca/micro_fr_pred/data/"


@TaskGenerator
def load_study(name: str, folder):
    return Study(name, path.join(folder, name))


@TaskGenerator
def reconstruct_mag(mag: str, mag_folder: str, reconstruction_folder: str):
    input = path.join(mag_folder, f"{mag}.fa.gz")
    output = path.join(reconstruction_folder, f"{mag}.xml")
    command = f"carve --dna {input} --output {output} -i M9 -g M9 -v"
    subprocess.check_call(command, shell=True)


study = load_study(study_name, data_folder)

for m in bvalue(study.mags):
    reconstruct_mag(
        m,
        path.join(data_folder, study_name, "mags"),
        os.makedirs(
            path.join(data_folder, study_name, "reconstructions"), exist_ok=True
        ),
    )
    print(f"Finished reconstruction for {m}")
