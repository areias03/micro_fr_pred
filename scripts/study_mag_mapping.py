import argparse
import os

import pandas as pd

SPIRE_PATH = "/work/microbiome/db/SPIRE/1/"
MAG_PATH = "/work/microbiome/db/SPIRE/1/spire_representative_genomes/"


def map_metadata_study(metadata: pd.DataFrame, mags: list, study: pd.DataFrame) -> list:
    # Match the metadata to the study samples
    study_meta = metadata[metadata["derived_from_sample"].isin(study["sample_id"])]
    # Filter out the MAGs that are not in the study metadata
    mag_names = [f.split(".fa.gz")[0] for f in mags]
    study_meta = study_meta[study_meta["genome_id"].isin(mag_names)]
    # Add the MAG filename to the metadata
    study_meta["mag_path"] = MAG_PATH + study_meta["genome_id"] + ".fa.gz"
    # Map the MAG to its corresponding sample in a new column
    study_meta["sample_mag_map"] = study_meta[
        ["mag_path", "derived_from_sample"]
    ].apply(lambda x: ",".join(x[x.notnull()]), axis=1)
    # Store list of target MAGs and mapping to reconstruct
    target_list = study_meta["sample_mag_map"].to_numpy()
    return target_list


def main():
    parser = argparse.ArgumentParser(
        description="Map target MAGs from study metadata to folder for reconstruction",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "input",
        metavar="INPUT",
        help="Input study metadata file for the mapping. Each study in SPIRE has a metadata file that can be used here.",
    )
    args = parser.parse_args()
    # Load the metadata
    metadata = pd.read_table(
        os.path.join("/work/microbiome/db/SPIRE/1/spire_v1_genome_metadata.tsv")
    )
    # Load the list of MAGs
    mags = [
        f for f in os.listdir(MAG_PATH) if os.path.isfile(os.path.join(MAG_PATH, f))
    ]
    input_df = pd.read_table(args.input)
    targets = map_metadata_study(metadata, mags, input_df)
    sum_targets = 0
    # Create file with list of target MAGs
    with open(
        "/work/microbiome/users/areiasca/functional_prediction/data/list_of_targets.txt",
        "w",
    ) as f:
        f.write("\n".join(targets))
        sum_targets += 1
    print(f"Number of targets: {len(targets)}")


if __name__ == "__main__":
    main()
