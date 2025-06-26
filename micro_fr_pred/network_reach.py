from typing import Dict, List, Union

import numpy as np
import polars as pl
from micom.workflows.results import GrowthResults


def network_reach(
    results: GrowthResults,
    taxa: Union[None, str, List[str]] = None,
) -> (Dict, Dict):
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
