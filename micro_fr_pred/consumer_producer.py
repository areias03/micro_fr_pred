import polars as pl
import numpy as np
from typing import Dict, Union
from micom.interaction import MES
from micom.workflows.results import GrowthResults


def _consumer_producer_score(reaction_scores, n_exchanges, total_exchanges, mode):
    if mode == "balanced":
        return sum(reaction_scores) * (n_exchanges / total_exchanges)
    elif mode == "unbalanced":
        return sum(reaction_scores)


def consumer_producer(
    results: GrowthResults,
    taxa: Union[None, str] = None,
    mode: str = "balanced",
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
        mode : str
            The balancing mode for the scores. 'balanced' if the scores are to be
    balanced by the ratio of exchanges or 'unblanced' if not

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
        i = 0
        for line in temp.iter_rows(named=True):
            temp[i, "MES"] = mes.filter(
                (mes["sample_id"] == line["sample_id"])
                & (mes["metabolite"] == line["metabolite"])
            )["MES"][0]
            i += 1
        temp = temp.with_columns(reaction_score=(pl.col("flux") * pl.col("MES")))
        scores[t] = _consumer_producer_score(
            temp["reaction_score"].to_list(), len(temp), len(exchanges), mode="balanced"
        )

    for i in scores.keys():
        if scores[i] > np.percentile(list(scores.values()), 75):
            classification[i] = "Producer"
        elif np.percentile(list(scores.values()), 25) > scores[i]:
            classification[i] = "Consumer"
        else:
            classification[i] = "Mixed"

    return scores, classification
