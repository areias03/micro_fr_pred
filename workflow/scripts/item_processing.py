from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import os.path as path
import polars as pl
import gzip
import glob.glob
from Bio import SeqIO


from spirepy import Study, Sample


def get_abundances(sample: Sample, folder: str):
    abundances = {}
    mag_folder = path.join(folder, "mags")
    depths = sample.get_contig_depths()
    for f in sorted(glob.glob(path.join(mag_folder, "*.fa.gz"))):
        with gzip.open(f, "rt") as handle:
            headers = set(rec.id for rec in SeqIO.parse(handle, "fasta"))
        abundance = depths.filter(pl.col("contigName").is_in(headers))[
            "totalAvgDepth"
        ].sum()
        abundances[f] = abundance
    return abundances


def generate_manifest(sample: Sample, folder: str):
    manif = []
    mag_folder = path.join(folder, "mags")
    reconstruction_folder = path.join(folder, "reconstructions")
    abun = get_abundances(sample)
    for genome in sample.get_mags().iter_rows(named=True):
        manif.append(
            [
                genome["spire_id"],
                genome["domain"],
                genome["phylum"],
                genome["class"],
                genome["order"],
                genome["family"],
                genome["genus"],
                genome["species"],
                path.join(reconstruction_folder, f"{genome['spire_id']}.xml"),
                genome["sample_id"],
                abun[path.join(mag_folder, f"{genome['spire_id']}.fa.gz")],
            ]
        )

    manifest = pl.DataFrame(
        manif,
        schema={
            "id": str,
            "kingdom": str,
            "phylum": str,
            "class": str,
            "order": str,
            "family": str,
            "genus": str,
            "species": str,
            "file": str,
            "sample_id": str,
            "abundance": float,
        },
        strict=False,
        orient="row",
    )
    manifest = manifest.with_columns(
        (pl.col("abundance") / pl.col("abundance").sum()).alias("abundance")
    )
    return manifest


def main(item, type, workers, out_folder):
    if type == "study":
        target = Study(item)
        results = []
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(generate_manifest, s, out_folder)
                for s in target.get_samples()
            ]
            for future in as_completed(futures):
                results.append(future.result())
        study_manifest = pl.concat(results)
        study_manifest.write_csv(path.join(out_folder, "manifest.csv"))

    elif type == "sample":
        target = Sample(item)
        sample_manifest = generate_manifest(target, out_folder)
        sample_manifest.write_csv(path.join(out_folder, "manifest.csv"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        metavar="INPUT",
        help="Input file (list of mags to download)",
        type=str,
    )
    parser.add_argument("-t", "--type", choices=["sample", "study"], default="study")
    parser.add_argument("-w", "--workers", dest="n_workers", type=int, default=1)
    parser.add_argument("-o", "--output", dest="output", help="output folder")
    args = parser.parse_args()
