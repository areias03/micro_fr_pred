from spirepy import Study
from spirepy.sample import Sample
import polars as pl


def sample_manifest(sample: Sample) -> pl.DataFrame:
    pass


def study_manifest(study: Study) -> pl.DataFrame:
    manifest = []
    for s in study.samples:
        manifest.append(sample_manifest(s))
    return pl.concat(manifest)
