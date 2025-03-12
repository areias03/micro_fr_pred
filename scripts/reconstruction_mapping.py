import os
from concurrent.futures import ThreadPoolExecutor

# import subprocess
import urllib.request
import io
import pandas as pd
import argparse
from micro_fr_pred.util import clean_emapper_data
import requests

data_folder = "/work/microbiome/users/areiasca/micro_fr_pred/data/"


def process_sample(samp):
    # Get metadata
    url = requests.get(f"https://spire.embl.de/api/sample/{samp}?format=tsv").text
    sample_meta = pd.read_csv(io.StringIO(url), sep="\t")
    print(
        f"Sample: {samp}\tNo. of MAGs: {len(sample_meta.spire_id)}\nList of MAGS:{sample_meta.spire_id.tolist()}"
    )
    # Download EggNOG-mapper data
    urllib.request.urlretrieve(
        f"https://spire.embl.de/download_eggnog/{samp}",
        f"{data_folder}{samp}/emapper_annotations.gz",
    )
    eggnog_data = clean_emapper_data(f"{data_folder}{samp}/emapper_annotations.gz")
    eggnog_data.to_csv(f"{data_folder}{samp}/emapper_annotations.tsv", sep="\t")
    # Map ,download and process MAGs
    mag_list = sample_meta.spire_id.tolist()
    for mag in mag_list:
        urllib.request.urlretrieve(
            f"https://spire.embl.de/download_file/{mag}",
            f"{data_folder}{samp}/mags/{mag}.fa.gz",
        )
        # Print reconstruction command
        command = f"carve --dna {data_folder}{samp}/mags/{mag}.fa.gz --egg {data_folder}{samp}/emapper_annotations.tsv -o {data_folder}{samp}/reconstructions/{mag}.xml"
        print(command)


def main():
    parser = argparse.ArgumentParser(
        description="Map target MAGs from study metadata to folder for reconstruction",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "study",
        help="Study ID from SPIRE to download data and reconstruct",
    )
    args = parser.parse_args()
    study = args.study
    study_folder = os.path.join(data_folder, study)
    os.makedirs(study_folder, exist_ok=True)
    stm = requests.get(f"https://spire.embl.de/api/study/{study}?format=tsv").text
    study_meta = pd.read_csv(io.StringIO(stm), sep="\t")
    list_of_samples = study_meta.sample_id.tolist()
    for sample in list_of_samples:
        sample_folder = f"{study_folder}/{sample}/"
        mag_folder = f"{sample_folder}/mags/"
        model_folder = f"{sample_folder}/reconstructions/"
        os.makedirs(sample_folder, exist_ok=True)
        os.makedirs(mag_folder, exist_ok=True)
        os.makedirs(model_folder, exist_ok=True)
    with ThreadPoolExecutor as executor:
        executor.map(process_sample, list_of_samples)


if __name__ == "__main__":
    main()
