from concurrent.futures import ThreadPoolExecutor
import polars as pl
import os
import os.path as path
import argparse
import urllib.request
import subprocess
import glob
import tempfile

def split_eggnog_data(mag, eggnog_data: str):
    with gzip.open(mag, "rt") as handle:
        headers = set(rec.id for rec in SeqIO.parse(handle, "fasta"))
    return eggnog_data.filter(pl.col("#query").is_in(headers))
       

def reconstruct_mag(mag, eggnog_data: str, reconstruction_folder: str):
    output = path.join(reconstruction_folder, f"{mag.split('/')[-1].strip('.fa.gz')}.xml")
    # command = f"carve --dna {mag} --output {output} -g M9 -v"
    command = f"carve --egg {eggnog_data} --output {output} -v"
    subprocess.check_call(command, shell=True)

def process_mag(mag, eggnog_data, reconstruction_folder):
    eggnog_path = path.join(reconstruction_folder, f"{mag}.csv")
    split_eggnog_data(mag, eggnog_data).write_csv(eggnog_path)
    reconstruct_mag(mag, eggnog_path, reconstruction_folder)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        metavar="INPUT",
        help="mags folder",
        type=str,
    )
    parser.add_argument("-o", "--output", dest="output", help="output folder")
    parser.add_argument("-e", "--egg", dest="eggnog", help="eggnog data")
    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    mags = sorted(glob.glob(path.join(args.input, "*.fa.gz")))
    egg = pl.read_csv(args.eggnog)
    with ThreadPoolExecutor() as executor:
        [executor.submit(process_mag, mag, egg, args.output) for mag in mags]
