from micro_fr_pred import Study
import os


data_folder = "/home/alexandre/Documents/projects/micro_fr_pred/data"
study_name = "Lloyd-Price_2019_HMP2IBD"

study = Study(study_name, os.path.join(data_folder, study_name))

# for sam in study.samples:
#     print(str(sam))


print(repr(study.samples[-1].metadata))
