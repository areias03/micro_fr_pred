import subprocess
from micro_fr_pred.logger import logger
import io
import pandas as pd
from micro_fr_pred.util import clean_emapper_data
import urllib.request
import os
import requests


class Study:
    """
    A study from SPIRE.

    This class represents a study from the SPIRE database. It automatically
    fetches metadata and automates the initialization of samples to further use
    to obtain its genomic, geographical or other types of data provided by it.

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
            logger.warning("No study metadata, downloading from SPIRE...\n")
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


class Sample:
    """
    A sample from SPIRE.

    This class represents a sample from the SPIRE database. It is designed to
    provide all the properties and methods to allow work with samples and
    provide tools for automation and scalability.

    Attributes:

    id: str
        Internal ID for each study.
    study: Study
        Study ID to which the sample belongs to.
    """

    def __init__(self, id: str, study: Study):
        """
        Creates a new sample object.
        """
        self.id = id
        self.study = study
        self.out_folder = f"{study.folder}{self.id}/"
        self._eggnog_data = None
        self._mags = None
        self._protein_bins = None
        self._metadata = None

        os.makedirs(self.out_folder, exist_ok=True)

    @property
    def eggnog_data(self):
        if self._eggnog_data is None:
            urllib.request.urlretrieve(
                f"https://spire.embl.de/download_eggnog/{self.id}",
                f"{self.out_folder}/emapper_annotations.gz",
            )
            eggnog_data = clean_emapper_data(
                f"{self.out_folder}/emapper_annotations.gz"
            )
            eggnog_data.to_csv(f"{self.out_folder}/emapper_annotations.tsv", sep="\t")
            self._eggnog_data = eggnog_data
        return self._eggnog_data

    @property
    def mags(self):
        if self._mags is None:
            self._mags = self.metadata.spire_id.tolist()
        return self._mags

    @property
    def protein_bins(self):
        return self._protein_bins

    @property
    def metadata(self):
        if self._metadata is None:
            url = requests.get(
                f"https://spire.embl.de/api/sample/{self.id}?format=tsv"
            ).text
            sample_meta = pd.read_csv(io.StringIO(url), sep="\t")
            self._metadata = sample_meta
        return self._metadata

    def download_mags(self):
        mag_folder = f"{self.out_folder}mags/"
        os.makedirs(mag_folder, exist_ok=True)
        for mag in self.mags:
            urllib.request.urlretrieve(
                f"https://spire.embl.de/download_file/{mag}",
                f"{mag_folder}{mag}.fa.gz",
            )

    def reconstruct(self):
        reconstruction_folder = f"{self.out_folder}reconstructions/"
        os.makedirs(reconstruction_folder, exist_ok=True)
        if len(os.listdir(f"{self.out_folder}mags/")) == 0:
            self.download_mags()
        command = f"carve -r {self.out_folder}mags/*.fa.gz --output {self.out_folder}reconstructions/"
        subprocess.check_call(command, shell=True)
