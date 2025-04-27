import os.path as path
from jug import TaskGenerator
from micro_fr_pred import Study, Sample

study_name = "Lloyd-Price_2019_HMP2IBD"
data_folder = "/work/microbiome/users/areiasca/micro_fr_pred/data/"


@TaskGenerator
def process_sample(sample: Sample):
    sample.download_mags()
    sample.reconstructions
    sample.manifest.to_csv(f"{sample.out_folder}sample_manifest.csv")


def main():
    study = Study(study_name, path.join(data_folder, study_name))
    for s in study.samples:
        process_sample(s)


if __name__ == "__main__":
    main()
