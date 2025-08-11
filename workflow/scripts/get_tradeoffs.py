import argparse
import os
from micom.workflows import build
from micom.workflows import tradeoff
from micom.qiime_formats import load_qiime_medium
from micom.viz import plot_tradeoff
import pandas as pd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        metavar="INPUT",
        help="Manifest file",
        type=str,
    )
    parser.add_argument("-m", "--models", dest="models", help="model folder")
    parser.add_argument(
        "-g", "--growth_medium", dest="medium", help="growth medium file"
    )
    parser.add_argument("-o", "--output", dest="output", help="output file prefix")
    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    data = pd.read_csv(args.input)
    manifest = build(
        data, model_db=None, out_folder=args.models, cutoff=0.0001, threads=24
    )
    medium = load_qiime_medium(args.medium)
    tradeoff_rates = tradeoff(
        manifest, model_folder=args.models, medium=medium, threads=24
    )
    tradeoff_rates.to_csv(os.path.join(args.output, "tradeoffs.csv"))
    plot_tradeoff(tradeoff_rates, filename=os.path.join(args.output, "tradeoffs.html"))
    tradeoff_rates.groupby("tradeoff").apply(
        lambda df: (df.growth_rate > 1e-6).sum()
    ).reset_index().to_csv(os.path.join(args.output, "tradeoffs_count.csv"))
