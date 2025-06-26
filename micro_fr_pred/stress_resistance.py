import random
from typing import Dict, List, Union

import numpy as np
import polars as pl
from micom.interaction import MES
from micom.workflows.results import GrowthResults


def stress_resistance(
    results: GrowthResults,
    taxa: Union[None, str, List[str]] = None,
) -> (Dict, Dict):
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
