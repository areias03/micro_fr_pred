from micom.qiime_formats import load_qiime_medium
from os.path import split, join


this_dir, _ = split(__file__)


def western_diet_gut():
    return load_qiime_medium(join(this_dir, "western_diet_gut.qza"))
