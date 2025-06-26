import random
from typing import Dict, List, Union

import numpy as np
import polars as pl
from micom.workflows.results import GrowthResults


def _metabolic_independence_score(taxon, samples, res):
    pass


def metabolic_independence(
    res: str,
    # results: GrowthResults,
    taxa: Union[None, str, List[str]] = None,
) -> (Dict, Dict):
    """Calculate the Metabolic Independece (MI) score for a taxon.

    This parameter is...

    Arguments
    ---------

    results : micom.workflows.results.GrowthResults
        The growth results to use.

    taxa : str, list of str, or None
        The focal taxa to use. Can be a single taxon, a list of taxa or None in which
        case all taxa are considered.

    Returns
    -------
    polars.DataFrame
        The scores for each taxon and their respectve classification.
    """
    growth_rates = pl.read_csv(res)
    # growth_rates = pl.from_pandas(results.growth_rates)

    taxon_list = growth_rates["taxon"].unique().to_list()

    for t in sorted(taxon_list):
        origin_sample = growth_rates.filter(growth_rates["taxon"] == t)[
            "sample_id"
        ].to_list()[0]
        random_samples = random.sample(
            sorted(
                [
                    s
                    for s in growth_rates["sample_id"].unique().to_list()
                    if s != origin_sample
                ]
            ),
            10,
        )
        print(
            f"Taxon: {t}\tSample: {origin_sample}\n Randomised samples:{random_samples}"
        )

    return None
