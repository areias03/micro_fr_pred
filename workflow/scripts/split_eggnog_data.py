from concurrent.futures import ThreadPoolExecutor
import os
import os.path as path
import argparse
import gzip
from Bio import SeqIO
import glob
import polars as pl
from util import get_ncpus


def split_eggnog_data(mag, eggnog_data, output_folder):
    with gzip.open(mag, "rt") as handle:
        headers = [rec.id for rec in SeqIO.parse(handle, "fasta")]
    print(mag, headers)
    split = eggnog_data.filter(
        (pl.col("#query").str.contains_any(headers))
        & (pl.col("sample") == mag.split("/")[-1].split("-")[0])
    )
    split.write_csv(
        path.join(
            output_folder, f"{mag.split('/')[-1].split('-')[1].strip('.fa.gz')}.tsv"
        ),
        separator="\t",
    )


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
    with ThreadPoolExecutor(max_workers=get_ncpus()) as executor:
        [executor.submit(split_eggnog_data, mag, egg, args.output) for mag in mags]
