import os
from concurrent.futures import ProcessPoolExecutor

import subprocess
import urllib.request
import io
import pandas as pd
import argparse
from micro_fr_pred.util import clean_emapper_data
from micro_fr_pred.ncpus import get_ncpus
import requests

data_folder = "/work/microbiome/users/areiasca/micro_fr_pred/data/"


def process_sample(args):
    samp = args[0]
    stu_fol = args[1]
    sample_folder = f"{stu_fol}/{samp}/"
    mag_folder = f"{sample_folder}/mags/"
    model_folder = f"{sample_folder}/reconstructions/"
    os.makedirs(sample_folder, exist_ok=True)
    os.makedirs(mag_folder, exist_ok=True)
    os.makedirs(model_folder, exist_ok=True)
    # Get metadata
    url = requests.get(f"https://spire.embl.de/api/sample/{samp}?format=tsv").text
    sample_meta = pd.read_csv(io.StringIO(url), sep="\t")
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
        command = f"carve --dna {sample_folder}mags/{mag}.fa.gz --output {sample_folder}reconstructions/{mag}.xml"
        subprocess.check_call(command, shell=True)


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
    with ProcessPoolExecutor(max_workers=get_ncpus()) as executor:
        arg1 = study_meta.sample_id.tolist()
        arg2 = list([str(study_folder)] * len(arg1))
        args = list(zip(arg1, arg2))
        tasks = executor.map(process_sample, args)
        # Iterating over the outputs will trigger error propagation
        # (If any of the tasks raised an exception, it will be re-raised here)
        for t in tasks:
            pass


if __name__ == "__main__":
    main()
