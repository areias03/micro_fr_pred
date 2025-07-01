import random
from typing import Dict, List, Union
import tempfile

import numpy as np
import polars as pl
from micom.workflows.results import GrowthResults
from micom.workflows import build, grow
from .data import western_diet_gut


def _randomize_samples(rates, origin, manif):
    random_samples = random.sample(
        sorted([s for s in rates["sample_id"].unique().to_list() if s != origin]),
        10,
    )
    random_manifest = manif.filter(manif["sample_id"].is_in(random_samples))
    return random_samples, random_manifest


def _compute_adjusted_abundances(sample, random_manif, abun):
    sample_manif = random_manif.filter(random_manif["sample_id"] == sample)
    adjusted_abun = abun / len(sample_manif)
    adjusted_manif = sample_manif.with_columns(pl.col("abundance") - adjusted_abun)
    return adjusted_manif


def _simulate_randomized_samples(manif):
    medium = western_diet_gut()
    with tempfile.TemporaryDirectory() as tmpdir:
        manifest = build(
            manif.to_pandas(),
            out_folder=tmpdir,
            model_db=None,
            cutoff=0.0001,
            threads=2,
        )
        res = grow(
            manifest, model_folder=tmpdir, medium=medium, tradeoff=0.5, threads=2
        )
        return res.growth_rates


def _metabolic_independence_score(original_growth_rate, randomized_growth_rates):
    return (original_growth_rate - np.mean(randomized_growth_rates)) / np.std(
        randomized_growth_rates
    )


def metabolic_independence(
    growth_rates: pl.DataFrame,
    manifest: pl.DataFrame,
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
    # growth_rates = pl.from_pandas(results.growth_rates)

    taxon_list = growth_rates["taxon"].unique().to_list()
    if taxa is not None:
        taxon_list = [t for t in taxon_list if t in taxa]

    res_list = {}
    for t in sorted(taxon_list):
        origin_sample = growth_rates.filter(growth_rates["taxon"] == t)[
            "sample_id"
        ].to_list()[0]
        abundance = growth_rates.filter(growth_rates["taxon"] == t)["abundance"][0]
        randomized_samples, randomized_manifest = _randomize_samples(
            growth_rates, origin_sample, manifest
        )
        manifest_list = []
        for s in randomized_samples:
            adjusted_manifest = _compute_adjusted_abundances(
                s, randomized_manifest, abundance
            )
            fitted_manifest = pl.concat(
                [adjusted_manifest, manifest.filter(manifest["id"] == t)]
            ).with_columns(pl.col("sample_id").str.replace(origin_sample, s))
            manifest_list.append(fitted_manifest)
        combined_manifest = pl.concat(manifest_list)
        simulated_rates = pl.from_pandas(
            _simulate_randomized_samples(combined_manifest)
        )
        randomized_rates = simulated_rates.filter(simulated_rates["taxon"] == t)[
            "growth_rate"
        ].to_list()
        original_rate = growth_rates.filter(growth_rates["taxon"] == t)["growth_rate"][
            0
        ]
        score = _metabolic_independence_score(original_rate, randomized_rates)
        print(
            f"Original: {original_rate}\nRandomized: {randomized_rates}\nScore: {_metabolic_independence_score(original_rate, randomized_rates)}"
        )
        res_list[t] = {
            "Original rate": original_rate,
            "Randomized rates": randomized_rates,
            "Score": score,
        }

    return res_list
