import os
import os.path as path
import urllib

from jug import TaskGenerator, bvalue

from spirepy import Study

study_name = "Lloyd-Price_2019_HMP2IBD"
data_folder = "/work/microbiome/users/areiasca/micro_fr_pred/data/"


@TaskGenerator
def load_study(name: str):
    return Study(name)


@TaskGenerator
def download_mag(mag: str, mag_folder: str):
    os.makedirs(mag_folder, exist_ok=True)
    print("Started downloading MAG:\t", mag)
    mag_out = path.join(mag_folder, f"{mag}.fa.gz")
    urllib.request.urlretrieve(
        f"https://spire.embl.de/download_file/{mag}",
        mag_out,
    )
    print(f"Finished downloading {mag}")
    return mag_out


study = load_study(study_name)

for m in bvalue(study).mags.iter_rows(named=True):
    download_mag(
        m["genome_id"],
        path.join(data_folder, study_name, m["derived_from_sample"], "mags"),
    )
