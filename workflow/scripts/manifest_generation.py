from concurrent.futures import ThreadPoolExecutor, as_completed
import argparse
import os.path as path
import polars as pl
import gzip
import glob
from Bio import SeqIO
from util import get_ncpus


from spirepy import Study, Sample


def get_abundances(sample: Sample, mag_folder):
    abundances = {}
    depths = sample.get_contig_depths()
    for f in sorted(
        [
            mag
            for mag in glob.glob(path.join(mag_folder, "*.fa.gz"))
            if mag.split("/")[-1].strip(".fa.gz")
            in sample.get_mags()["genome_id"].to_list()
        ]
    ):
        with gzip.open(f, "rt") as handle:
            headers = set(rec.id for rec in SeqIO.parse(handle, "fasta"))
        abundance = depths.filter(pl.col("contigName").is_in(headers))[
            "totalAvgDepth"
        ].sum()
        abundances[f] = abundance
    return abundances


def generate_manifest(sample: Sample, mag_folder: str, reconstruction_folder: str):
    manif = []
    abun = get_abundances(sample, mag_folder)
    for genome in sample.get_mags().iter_rows(named=True):
        manif.append(
            [
                genome["genome_id"],
                genome["domain"],
                genome["phylum"],
                genome["class"],
                genome["order"],
                genome["family"],
                genome["genus"],
                genome["species"],
                path.join(reconstruction_folder, f"{genome['genome_id']}.xml"),
                genome["derived_from_sample"],
                abun[path.join(mag_folder, f"{genome['genome_id']}.fa.gz")],
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


def main(item, type, out_file, mag_folder, reconstruction_folder):
    if type == "study":
        target = Study(item)
        results = []
        with ThreadPoolExecutor(max_workers=get_ncpus()) as executor:
            futures = [
                executor.submit(generate_manifest, s, mag_folder, reconstruction_folder)
                for s in target.get_samples()
            ]
            for future in as_completed(futures):
                results.append(future.result())
        study_manifest = pl.concat(results)
        study_manifest.write_csv(out_file)

    elif type == "sample":
        target = Sample(item)
        sample_manifest = generate_manifest(target, mag_folder, reconstruction_folder)
        sample_manifest.write_csv(out_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "item",
        metavar="ITEM",
        help="Input item (study or sample ID)",
        type=str,
    )
    parser.add_argument("-t", "--type", choices=["sample", "study"], default="study")
    parser.add_argument("-o", "--output", dest="output", help="output folder")
    parser.add_argument("-m", "--mags", dest="mags", help="mags folder")
    parser.add_argument(
        "-r", "--reconstruction", dest="reconstructions", help="reconstructions folder"
    )
    args = parser.parse_args()
    main(args.item, args.type, args.output, args.mags, args.reconstructions)
