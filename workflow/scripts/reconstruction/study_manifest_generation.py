import glob
import gzip
import os.path as path

import polars as pl
from Bio import SeqIO
from spirepy import Study
from spirepy.sample import Sample
from tqdm import tqdm

study_name = "Lloyd-Price_2019_HMP2IBD"
data_folder = "/home/alexandre/Documents/projects/micro_fr_pred/data/"


def get_abundances(sample: Sample):
    abundances = {}
    mag_folder = path.join(data_folder, study_name, sample.id, "mags")
    depths = pl.read_csv(
        # f"/work/microbiome/global_data_spire/SPIRE/studies/{study_name}/psa_megahit/psb_metabat2/{sample.id}_aligned_to_{sample.id}.depths",
        f"/home/alexandre/microbiome/global_data_spire/SPIRE/studies/{study_name}/psa_megahit/psb_metabat2/{sample.id}_aligned_to_{sample.id}.depths",
        separator="\t",
    )
    for f in sorted(glob.glob(path.join(mag_folder, "*.fa.gz"))):
        with gzip.open(f, "rt") as handle:
            headers = set(rec.id for rec in SeqIO.parse(handle, "fasta"))
        abundance = depths.filter(depths["contigName"].is_in(headers))[
            "totalAvgDepth"
        ].sum()
        abundances[f] = abundance
    return abundances


def generate_manifest(sample: Sample):
    manif = []
    mag_folder = path.join(data_folder, study_name, sample.id, "mags")
    reconstruction_folder = path.join(
        data_folder, study_name, sample.id, "reconstructions"
    )
    abun = get_abundances(sample)
    for genome in sample.mags.iter_rows(named=True):
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
    oname = path.join(data_folder, study_name, sample.id, "sample_manifest.csv")
    manifest.write_csv(oname)
    return manifest


def load_study(study: str):
    return Study(study)


def merge_partials(partials):
    study_manifest = pl.concat(partials, how="vertical")
    study_manifest = study_manifest.group_by("sample_id")
    oname = path.join(data_folder, study_name, "study_manifest.csv")
    study_manifest.write_csv(oname)
    return oname


def main():
    study = load_study(study_name)

    manifest_list = []
    for s in tqdm(study.samples):
        print(s.id)
        try:
            manifest = generate_manifest(s)
            manifest_list.append(manifest)
        except:
            print("Passed")
            pass
    merge_partials(manifest_list)


if __name__ == "__main__":
    main()
