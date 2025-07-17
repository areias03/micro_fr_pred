from spirepy.sample import Sample
import os.path as path
import gzip
import polars as pl
import micom
from Bio import SeqIO
from glob import glob

study_name = "Lloyd-Price_2019_HMP2IBD"
data_folder = "/work/microbiome/users/areiasca/micro_fr_pred/data/"

def get_abundances(sample: Sample):
    abundances = {}
    mag_folder = path.join(data_folder, study_name, sample.id, "mags")
    depths = pl.read_csv(
        f"/work/microbiome/global_data_spire/SPIRE/studies/{study_name}/psa_megahit/psb_metabat2/{sample.id}_aligned_to_{sample.id}.depths",
        separator="\t",
    )
    for f in sorted(glob(path.join(mag_folder, "*.fa.gz"))):
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
    return manifest

def 
