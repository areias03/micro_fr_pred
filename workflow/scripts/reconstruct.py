import argparse
import glob
import os
import os.path as path
import subprocess
from concurrent.futures import ThreadPoolExecutor
from util import get_ncpus


def reconstruct_mag(eggnog_data, reconstruction_folder: str):
    output = path.join(
        reconstruction_folder, f"{eggnog_data.split('/')[-1].removesuffix('.tsv')}.xml"
    )
    # command = f"carve --dna {mag} --output {output} -g M9 -v"
    command = f"carve --egg {eggnog_data} --output {output} -v"
    subprocess.check_call(command, shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        metavar="INPUT",
        help="eggnog data folder",
        type=str,
    )
    parser.add_argument("-o", "--output", dest="output", help="output folder")
    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    eggnog = sorted(glob.glob(path.join(args.input, "*.tsv")))
    print(eggnog)
    with ThreadPoolExecutor(max_workers=get_ncpus()) as executor:
        [executor.submit(reconstruct_mag, egg, args.output) for egg in eggnog]
