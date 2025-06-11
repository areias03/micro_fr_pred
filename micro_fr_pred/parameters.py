import random
from typing import Dict, List, Union

import numpy as np
import polars as pl
from micom.interaction import MES
from micom.workflows.results import GrowthResults


def _consumer_producer_score(taxon):
    pass


def consumer_producer(
    results: GrowthResults, taxa: Union[None, str] = None
) -> (Dict, Dict):
    """Calculate the consumer/producer score for a taxon.

    This parameter is defined as the sum of all import and export fluxes
    multiplied by their repective metabolite's Metabolite Exchange Score (MES).
    It represents the harmonic mean of all fluxes for a single taxon. A negative
    value indicates that the taxon consumes more impactful metabolites than what
    it produces. A positive value indicates higher rates of production of impactful
    metabolites to the whole community.

    Arguments
    ---------

    results : micom.workflows.results.GrowthResults
        The growth results to use.

    taxa : str or None
        The focal taxa to use. Can be a single taxon or None in which
        case all taxa are considered.

    Returns
    -------
    polars.DataFrame
        The scores for each taxon.

    polars.DataFrame
        The classification for each taxon.

    References
    ----------
    .. [1] Marcelino, V.R., et al.
           Disease-specific loss of microbial cross-feeding interactions in the human gut
           Nat Commun 14, 6546 (2023). https://doi.org/10.1038/s41467-023-42112-w


    """
    exchanges = pl.from_pandas(results.exchanges)
    mes = pl.from_pandas(MES(results))
    exchanges = exchanges.filter(exchanges["taxon"] != "medium")

    if taxa is not None:
        exchanges = exchanges.filter(exchanges["taxon"] == taxa)
    exchanges = exchanges.with_columns(MES=pl.lit(0))

    taxon_list = exchanges["taxon"].unique().to_list()

    scores = {}
    classification = {}

    for t in sorted(taxon_list):
        temp = exchanges.filter(exchanges["taxon"] == t)
        size = len(temp)
        i = 0
        for line in temp.iter_rows(named=True):
            temp[i, "MES"] = mes.filter(
                (mes["sample_id"] == line["sample_id"])
                & (mes["metabolite"] == line["metabolite"])
            )["MES"][0]
            i += 1
        temp = temp.with_columns(
            score=(pl.col("flux") * pl.col("MES")) * (size / len(exchanges))
        )
        temp = temp.group_by("taxon").sum()

        scores[t] = temp["score"][0]

    for i in scores.keys():
        if scores[i] > np.percentile(list(scores.values()), 75):
            classification[i] = "Producer"
        elif np.percentile(list(scores.values()), 25) > scores[i]:
            classification[i] = "Consumer"
        else:
            classification[i] = "Mixed"

    return scores, classification


def metabolic_independence(
    res: str,
    # results: GrowthResults,
    taxa: Union[None, str, List[str]] = None,
) -> pl.DataFrame:
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

    References
    ----------
    .. [1] Watson, A.R., Füssel, J., Veseli, I. et al.
           Metabolic independence drives gut microbial colonization and resilience in health and disease.
           Genome Biol 24, 78 (2023). https://doi.org/10.1186/s13059-023-02924-x
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
                    for s in growth_rates.filter(growth_rates["taxon"] == t)
                    if s != origin_sample
                ]
            ),
            10,
        )
        print(
            f"Taxon: {t}\tSample: {origin_sample}\n Randomised samples:{random_samples}\tMistake in sampling: {origin_sample in random_samples}"
        )

    return None


def community_importance(
    results: GrowthResults,
    taxa: Union[None, str, List[str]] = None,
) -> pl.DataFrame:
    """Calculate and classify the impact a single microbe has on the community.

    This parameter

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

    References
    ----------
    .. [1] Marcelino, V.R., et al.
           Disease-specific loss of microbial cross-feeding interactions in the human gut
           Nat Commun 14, 6546 (2023). https://doi.org/10.1038/s41467-023-42112-w


    """
    return None


def network_reach(
    results: GrowthResults,
    taxa: Union[None, str, List[str]] = None,
) -> pl.DataFrame:
    """Calculate the reach of a microbe's metabolic network.

    This parameter

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
    return None


def stress_resistance(
    results: GrowthResults,
    taxa: Union[None, str, List[str]] = None,
) -> pl.DataFrame:
    """Calculate the resistance to metabolic stress.

    This parameter

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

    References
    ----------
    .. [1] Marcelino, V.R., et al.
           Disease-specific loss of microbial cross-feeding interactions in the human gut
           Nat Commun 14, 6546 (2023). https://doi.org/10.1038/s41467-023-42112-w


    """
    return None
