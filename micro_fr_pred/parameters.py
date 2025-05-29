import polars as pl
from micom.workflows.results import GrowthResults
from micom.interaction import MES
from typing import Union, List


def metabolic_independence(
    results: GrowthResults,
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
    return None


def consumer_producer(
    results: GrowthResults,
    taxa: Union[None, str, List[str]] = None,
) -> pl.DataFrame:
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
    exchanges = pl.from_pandas(results.exchanges)
    rates = pl.from_pandas(results.growth_rates)
    annotations = pl.from_pandas(results.annotations)
    mes = MES(results)
    mes = pl.from_pandas(mes)

    if taxa is not None:
        if taxa is str:
            exchanges = exchanges.filter(exchanges["taxon"] == taxa)
        else:
            exchanges = exchanges.filter(exchanges["taxon"].is_in(taxa))
    exchanges = exchanges.with_columns(MES=pl.lit(0))
    exchanges = exchanges.filter(exchanges["taxon"] != "medium")
    i = 0
    for line in exchanges.iter_rows(named=True):
        exchanges[i, "MES"] = mes.filter(
            (mes["sample_id"] == line["sample_id"])
            & (mes["metabolite"] == line["metabolite"])
        )["MES"][0]
        i += 1

    return exchanges


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
