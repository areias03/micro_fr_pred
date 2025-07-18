from concurrent.futures import ThreadPoolExecutor
import os
import os.path as path
import argparse
import urllib.request
from util import get_ncpus


def download_mag(mag, folder):
    mag_out = path.join(folder, f"{mag}.fa.gz")
    urllib.request.urlretrieve(
        f"https://spire.embl.de/download_file/{mag}",
        mag_out,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        metavar="INPUT",
        help="Input file (list of mags to download)",
        type=argparse.FileType("r"),
    )
    parser.add_argument("-o", "--output", dest="output", help="output folder")
    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    with args.input as f:
        mags = f.read().splitlines()
        with ThreadPoolExecutor(max_workers=get_ncpus()) as executor:
            [executor.submit(download_mag, mag, args.output) for mag in mags]
