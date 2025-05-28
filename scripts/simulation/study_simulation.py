import os.path as path

import pandas as pd

# from micom.data import test_medium
from micom.interaction import interactions, summarize_interactions
from micom.qiime_formats import load_qiime_medium
from micom.viz import (
    plot_exchanges_per_taxon,
    plot_focal_interactions,
    plot_growth,
    plot_mes,
)
from micom.workflows import build, grow
# from micom.workflows import complete_community_medium

study_name = "Lloyd-Price_2019_HMP2IBD"
# data_folder = "/home/alexandre/Documents/projects/micro_fr_pred/data/"
data_folder = "/work/microbiome/users/areiasca/micro_fr_pred/data/"


def main():
    out_folder = path.join(data_folder, study_name, "simulations")
    tax = pd.read_csv(path.join(data_folder, study_name, "study_manifest.csv"))
    manifest = build(
        taxonomy=tax,
        out_folder=out_folder,
        model_db=None,
        cutoff=0.0001,
        threads=32,
    )
    # Western diet medium
    medium = load_qiime_medium(path.join(data_folder, "western_diet_gut.qza"))
    res = grow(
        manifest,
        model_folder=out_folder,
        medium=medium,
        tradeoff=1.0,
        threads=32,
    )
    res.save(path.join(out_folder, "simulation_results.zip"))
    plot = plot_growth(res, filename=path.join(out_folder, "growth_rates.html"))
    plot = plot_mes(res, filename=path.join(out_folder, "mes.html"))
    plot = plot_exchanges_per_taxon(res, filename=path.join(out_folder, "niche.html"))
    full = interactions(res, taxa=None, threads=32)
    full.to_csv(path.join(out_folder, "interactions.csv"))
    plot = plot_focal_interactions(
        res, taxon=None, filename=path.join(out_folder, "niche.html")
    )
    summary = summarize_interactions(full)
    summary.to_csv(path.join(out_folder, "interactions_summary.csv"))


if __name__ == "__main__":
    main()
