import os
from concurrent.futures import ProcessPoolExecutor

import subprocess
import urllib.request
import io
import pandas as pd
import argparse
from micro_fr_pred.util import clean_emapper_data
import requests

data_folder = "./data/"


def process_sample(samp, stu_fol):
    sample_folder = f"{stu_fol}/{samp}/"
    mag_folder = f"{sample_folder}/mags/"
    model_folder = f"{sample_folder}/reconstructions/"
    os.makedirs(sample_folder, exist_ok=True)
    os.makedirs(mag_folder, exist_ok=True)
    os.makedirs(model_folder, exist_ok=True)
    # Get metadata
    url = requests.get(f"https://spire.embl.de/api/sample/{samp}?format=tsv").text
    sample_meta = pd.read_csv(io.StringIO(url), sep="\t")
    print(
        f"Sample: {samp}\tNo. of MAGs: {len(sample_meta.spire_id)}\nList of MAGS:{sample_meta.spire_id.tolist()}"
    )
    # Download EggNOG-mapper data
    urllib.request.urlretrieve(
        f"https://spire.embl.de/download_eggnog/{samp}",
        f"{sample_folder}emapper_annotations.gz",
    )
    eggnog_data = clean_emapper_data(f"{sample_folder}emapper_annotations.gz")
    eggnog_data.to_csv(f"{sample_folder}emapper_annotations.tsv", sep="\t")
    # Map ,download and process MAGs
    mag_list = sample_meta.spire_id.tolist()
    for mag in mag_list:
        urllib.request.urlretrieve(
            f"https://spire.embl.de/download_file/{mag}",
            f"{sample_folder}mags/{mag}.fa.gz",
        )
        # Print reconstruction command
        command = f"carve --dna {sample_folder}mags/{mag}.fa.gz --egg {sample_folder}emapper_annotations.tsv -o {sample_folder}reconstructions/{mag}.xml"
        print(command)
        subprocess.check_call(command)


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
    with ProcessPoolExecutor() as executor:
        tasks = executor.map(process_sample, study_meta.sample_id.tolist())
        # Iterating over the outputs will trigger error propagation
        # (If any of the tasks raised an exception, it will be re-raised here)
        for t in tasks:
            pass


if __name__ == "__main__":
    main()
