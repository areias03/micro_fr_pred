import argparse
import glob
import os
import os.path as path
import subprocess
from concurrent.futures import ThreadPoolExecutor


def reconstruct_mag(mag, eggnog_folder, reconstruction_folder: str):
    eggnog_data = path.join(eggnog_folder, f"{mag.split('/')[-1].strip('.fa.gz')}.tsv")
    output = path.join(
        reconstruction_folder, f"{mag.split('/')[-1].strip('.fa.gz')}.xml"
    )
    # command = f"carve --dna {mag} --output {output} -g M9 -v"
    command = f"carve --egg {eggnog_data} --output {output} -v"
    subprocess.check_call(command, shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        metavar="INPUT",
        help="mags folder",
        type=str,
    )
    parser.add_argument("-o", "--output", dest="output", help="output folder")
    parser.add_argument("-e", "--egg", dest="eggnog", help="eggnog data folder")
    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    mags = sorted(glob.glob(path.join(args.input, "*.fa.gz")))
    with ThreadPoolExecutor() as executor:
        [
            executor.submit(reconstruct_mag, mag, args.eggnog, args.output)
            for mag in mags
        ]
