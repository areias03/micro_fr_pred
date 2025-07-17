from concurrent.futures import ThreadPoolExecutor
import os
import os.path as path
import argparse
import gzip
from Bio import SeqIO
import glob
import polars as pl


def split_eggnog_data(mag, eggnog_data, output_folder):
    with gzip.open(mag, "rt") as handle:
        headers = [f"{rec.id}" for rec in SeqIO.parse(handle, "fasta")]
        for h in headers:
            headers.append(f"{h}_1")
            headers.append(f"{h}_2")
    print(mag, headers)
    print(eggnog_data.filter(eggnog_data["#query"].is_in(headers)))
    eggnog_data.filter(pl.col("#query").is_in(headers)).write_csv(
        path.join(output_folder, f"{mag.split('/')[-1].strip('.fa.gz')}.tsv"),
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
    print(egg.group_by("#query").len().describe())
    with ThreadPoolExecutor() as executor:
        [executor.submit(split_eggnog_data, mag, egg, args.output) for mag in mags]
