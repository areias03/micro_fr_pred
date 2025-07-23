import argparse
from micom.workflows import build
from micom.workflows import tradeoff
from micom.qiime_formats import load_qiime_medium
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        metavar="INPUT",
        help="Manifest file",
        type=str,
    )
    parser.add_argument("-o", "--output", dest="output", help="output file prefix")
    args = parser.parse_args()
    data = pd.read_csv(args.input)
    manifest = build(data, out_folder="models", cutoff=0.0001, threads=32)
    medium = load_qiime_medium("../data/western_diet_gut.qza")
    tradeoff_rates = tradeoff(
        manifest, model_folder="models", medium=medium, threads=32
    )
    tradeoff_rates.groupby("tradeoff").apply(
        lambda df: (df.growth_rate > 1e-6).sum()
    ).reset_index().to_csv(args.ouput)
