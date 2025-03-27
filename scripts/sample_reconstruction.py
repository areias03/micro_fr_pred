from micro_fr_pred import Study
import os


data_folder = "/home/alexandre/Documents/projects/micro_fr_pred/data"
study_name = "Lloyd-Price_2019_HMP2IBD"

study = Study(study_name, os.path.join(data_folder, study_name))


# Tests

print("###Study ###")
print(f"Name: {study.name}\n")
print(f"Metadata:\n {study.metadata}\n")

print("###Sample ###")
print(f"Name: {study.samples[-1].id}\n")
print(f"MAGs:\n {study.samples[-1].mags}\n")
print(f"Manifest:\n {study.samples[-1].manifest}\n")
