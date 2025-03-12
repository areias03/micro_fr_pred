import os
import io
import pandas as pd
import argparse
# from micro_fr_pred.util import clean_emapper_data

# import numpy as np
# from Bio import SeqIO
import requests

eggnog_url = "https://spire.embl.de/download_eggnog/"
protein_assembly_url = "https://swifter.embl.de/~fullam/spire/genes_per_study/Lloyd-Price_2019_HMP2IBD_spire_v1_genecalls_faa.tar"
mag_url = "https://swifter.embl.de/~fullam/spire/compiled/Lloyd-Price_2019_HMP2IBD_spire_v1_MAGs.tar"
study_metadata = "https://spire.embl.de/api/study/Lloyd-Price_2019_HMP2IBD?format=tsv"
sample_metadata = "https://spire.embl.de/api/sample/SAMN07510031?format=tsv"
data_folder = "../data/"


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
        os.makedirs(sample_folder, exist_ok=True)
        sam = requests.get(f"https://spire.embl.de/api/sample/{sample}?format=tsv").text
        sample_meta = pd.read_csv(io.StringIO(sam), sep="\t")
        print(f"Sample: {sample}\tNo. of MAGs: {len(sample_meta.spire_id)}")


if __name__ == "__main__":
    main()
