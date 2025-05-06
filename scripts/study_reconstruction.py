import os
import polars as pl
import os.path as path
import subprocess

from jug import TaskGenerator, bvalue, iteratetask

from spirepy import Study

study_name = "Lloyd-Price_2019_HMP2IBD"
# data_folder = "/work/microbiome/users/areiasca/micro_fr_pred/data/"
data_folder = "/home/alexandre/Documents/projects/micro_fr_pred/data/"


@TaskGenerator
def load_study(name: str):
    return Study(name)


@TaskGenerator
def generate_mag_list(study: Study):
    genome_metadata = pl.read_csv(
        path.join(data_folder, "spire_v1_genome_metadata.tsv.gz"), separator="\t"
    )
    return genome_metadata.filter(
        genome_metadata["genome_id"].is_in(
            bvalue(study).metadata["sample_id"].to_list()
        )
    )["genome_id"].to_list()


@TaskGenerator
def reconstruct_mag(mag: str, mag_folder: str, reconstruction_folder: str):
    print(f"Starting reconstruction for {m}")
    input = path.join(mag_folder, f"{mag}.fa.gz")
    output = path.join(reconstruction_folder, f"{mag}.xml")
    command = f"carve --dna {input} --output {output} -i M9 -g M9 -v"
    subprocess.check_call(command, shell=True)
    print(f"Finished reconstruction for {m}")


study = load_study(study_name)
mags = generate_mag_list(study)
print(mags)

for m in bvalue(mags):
    reconstruct_mag(
        m,
        path.join(data_folder, study_name, "mags"),
        os.makedirs(
            path.join(data_folder, study_name, "reconstructions"), exist_ok=True
        ),
    )
