import pandas as pd
import os
import numpy as np


def clean_emapper_data(file: str) -> pd.DataFrame:
    data = pd.read_csv(file, skipfooter=3, skiprows=4, sep="\t")
    data.columns = data.columns.str.replace("#", "")
    return data
