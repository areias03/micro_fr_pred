from study import Study
import os
import requests


class Sample(object):
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

    def __init__(self, id, out_folder=None):
        """
        Creates a new sample object.
        """
        self.id = id
        self.out_folder = out_folder
        self.eggnog_data = None
        self.mags = None
        self.protein_bins = None
        self.metadata = None
        self.study: Study = None

        if out_folder:
            if os.path.exists(out_folder):
                pass
            else:
                os.makedirs(out_folder)
        else:
            out_folder = f"./{id}"
            os.makedirs(out_folder)

    @property
    def protein_bins(self):
        return self.protein_bins

    @property
    def eggnog_data(self):
        return self.eggnog_data

    @property
    def genomic_asembly(self):
        return self.genomic_asembly
