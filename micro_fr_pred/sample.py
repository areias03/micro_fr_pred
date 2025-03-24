from study import Study
import subprocess
import io
import pandas as pd
from util import clean_emapper_data
import urllib.request
import os
import requests


class Sample:
    """
    A sample from SPIRE.

    This class represents a class from the SPIRE database. It is designed to
    provide all the propertuies and methods that allow to work with the samples and
    provide tools for automation and scalability.

    Attributes:

    id: str
        Internal ID for each study.
    out_folder: str, optional
        Output folder where everything will be downloaded and
        reconstructed. Will default to new folder in current directory with the name of
        the study.
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
