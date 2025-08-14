import argparse
import os
from micom.workflows import build, grow
from micom.interaction import interactions, MES
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
    parser.add_argument("-m", "--models", dest="models", help="model folder")
    parser.add_argument(
        "-g", "--growth_medium", dest="medium", help="growth medium file"
    )
    parser.add_argument("-t", "--tradeoff", dest="tradeoff", help="tradeoff value")
    parser.add_argument("-o", "--output", dest="output", help="output folder")
    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    data = pd.read_csv(args.input)
    manifest = build(
        data,
        model_db=None,
        out_folder=args.models,
        cutoff=0.0001,
        threads=4,
    )
    medium = load_qiime_medium(args.medium)
    res = grow(
        manifest,
        model_folder=args.models,
        medium=medium,
        tradeoff=args.tradeoff,
        threads=4,
    )
    res.growth_rates.to_csv(os.path.join(args.output, "growth_rates.csv"))
    res.exchanges.to_csv(os.path.join(args.output, "exchanges.csv"))
    res.annotations.to_csv(os.path.join(args.output, "annotations.csv"))
    interactions(res, taxa=None, threads=24).to_csv(
        os.path.join(args.output, "interactions.csv")
    )
    MES(res).to_csv(os.path.join(args.output, "mes.csv"))
