import io
import os
import subprocess
import urllib.request

import pandas as pd
import requests

from micro_fr_pred.logger import logger
from micro_fr_pred.util import clean_emapper_data
import gzinfo

# TODO: Add higher level class for overall SPIRE to help deal with queries and dealing with metadata.


def genome_metadata():
    data = pd.read_csv("data/spire_v1_genome_metadata.tsv.gz", sep="\t")
    return data


class Study:
    """
    A study from SPIRE.

    This class represents a study from the SPIRE database. It automatically
    fetches metadata and automates the initialization of samples to further use
    to obtain its genomic, geographical or other types of data provided by it.

    Attributes:

    name: str
        Internal ID for the study.
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
            sample_list = []
            for s in self.metadata.sample_id.tolist():
                sample = Sample(s, self)
                sample_list.append(sample)
            self._samples = sample_list
        return self._samples


class Sample:
    """
    A sample from SPIRE.

    This class represents a sample from the SPIRE database. It is designed to
    provide all the properties and methods to allow work with samples and
    provide tools for automation and scalability.

    Attributes:

    id: str
        Internal ID for the sample.
    study: Study
        Study ID to which the sample belongs to.
    """

    def __init__(self, id: str, study: Study):
        """
        Creates a new sample object.
        """
        self.id = id
        self.study = study
        self.out_folder = f"{study.folder}/{self.id}/"
        self._eggnog_data = None
        self._mags = None
        self._metadata = None
        self._manifest = None

        os.makedirs(self.out_folder, exist_ok=True)

    def __str__(self):
        return (
            f"This sample's id is {self.id} and belongs to the {self.study.name} study."
        )

    def __repr__(self):
        return f"Sample id: {self.id} \tStudy: {self.study.name}"

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
    def mags(self, download: bool = False):
        if self._mags is None:
            spire_meta = genome_metadata()
            masked = spire_meta.loc[spire_meta["derived_from_sample"] == self.id]
            self._mags = masked
        if download:
            self.download_mags()
        return self._mags

    @property
    def metadata(self):
        if self._metadata is None:
            logger.warning("No sample metadata, downloading from SPIRE...\n")
            url = requests.get(
                f"https://spire.embl.de/api/sample/{self.id}?format=tsv"
            ).text
            sample_meta = pd.read_csv(io.StringIO(url), sep="\t")
            self._metadata = sample_meta
        return self._metadata

    @property
    def manifest(self):
        if self._manifest is None:
            self._manifest = self.generate_manifest()
        return self._manifest

    def download_mags(self):
        mag_folder = f"{self.out_folder}mags/"
        os.makedirs(mag_folder, exist_ok=True)
        for mag in self.mags:
            urllib.request.urlretrieve(
                f"https://spire.embl.de/download_file/{mag}",
                f"{mag_folder}{mag}.fa.gz",
            )

        for f in os.listdir(mag_folder):
            filename = gzinfo.read_gz_info(os.path.join(mag_folder, f))
            print(filename.fname)

    def generate_manifest(self):
        manif = []
        reconstruction_folder = f"{self.out_folder}reconstructions/"
        for _, genome in self.mags.iterrows():
            manif.append(
                [
                    genome.genome_id,
                    genome.genus,
                    genome.species,
                    f"{reconstruction_folder}{genome.genome_id}.xml",
                    genome.derived_from_sample,
                    0,
                ]
            )

        manifest = pd.DataFrame(
            manif, columns=["id", "genus", "species", "file", "sample_id", "abundance"]
        )
        manifest.groupby("sample_id")
        return manifest

    def reconstruct(self):
        reconstruction_folder = f"{self.out_folder}reconstructions/"
        os.makedirs(reconstruction_folder, exist_ok=True)
        if (
            not os.path.exists(f"{self.out_folder}mags/")
            or len(os.listdir(f"{self.out_folder}mags/")) == 0
        ):
            self.download_mags()
        command = f"carve -r {self.out_folder}mags/*.fa.gz --output {self.out_folder}reconstructions/"
        subprocess.check_call(command, shell=True)
