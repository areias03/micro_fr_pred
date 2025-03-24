import os
import io
import pandas as pd
import requests
from sample import Sample


class Study:
    """
    A study from SPIRE.

    This class represents...

    Attributes:

    name: str
        And description.
    """

    def __init__(self, name: str, out_folder: str):
        self.name = name
        # self.sample = list
        self.folder = out_folder
        self._metadata = None
        self._samples = None

        os.makedirs(self.folder, exist_ok=True)

    @property
    def metadata(self):
        if self._metadata is None:
            print("No metadata, downloading from SPIRE")
            url = requests.get(
                f"https://spire.embl.de/api/study/{self.name}?format=tsv"
            ).text
            study_meta = pd.read_csv(io.StringIO(url), sep="\t")
            self._metadata = study_meta
        return self._metadata

    @property
    def samples(self):
        if self._samples is None:
            self._samples = self.metadata.sample_id.tolist()
        return self._samples
