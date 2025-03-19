import micom
# import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import io
import requests

genome_metadata = pd.read_csv(
    "~/microbiome/users/areiasca/micro_fr_pred/data/spire_v1_genome_metadata.tsv.gz",
    sep="\t",
)

study = "Lloyd-Price_2019_HMP2IBD"
stm = requests.get(f"https://spire.embl.de/api/study/{study}?format=tsv").text
study_meta = pd.read_csv(io.StringIO(stm), sep="\t")
mask = [
    sample in study_meta.sample_id.tolist()
    for sample in genome_metadata["derived_from_sample"]
]


def main():
    manif = []
    unmapped = 0
    mapped = 0
    for _, genome in genome_metadata[mask].iterrows():
        manif.append(
            [
                genome.genome_id,
                genome.genus,
                genome.species,
                f"/work/microbiome/users/areiasca/micro_fr_pred/data/{study}/{genome.derived_from_sample}/reconstructions/{genome.genome_id}.xml",
                genome.derived_from_sample,
                0,
            ]
        )
        if pd.isna(genome.species):
            unmapped += 1
        else:
            mapped += 1
    manifest = pd.DataFrame(
        manif, columns=["id", "genus", "species", "file", "sample_id", "abundance"]
    )
    manifest.groupby("sample_id")
    print(
        f"\nUnmapped: {unmapped}\{unmapped + mapped}\tPercentage: {round(((unmapped / mapped) * 100), 2)}%"
    )


if __name__ == "__main__":
    main()
