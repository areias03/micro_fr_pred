from concurrent.futures import ThreadPoolExecutor
import os
import os.path as path
import argparse
import urllib.request
import subprocess


def download_mag(mag, folder):
    mag_out = path.join(folder, f"{mag}.fa.gz")
    urllib.request.urlretrieve(
        f"https://spire.embl.de/download_file/{mag}",
        mag_out,
    )


def reconstruct_mag(mag: str, mag_folder: str, reconstruction_folder: str):
    os.makedirs(reconstruction_folder, exist_ok=True)
    print(f"Started reconstructing:\t {mag}")
    input = path.join(mag_folder, f"{mag}.fa.gz")
    output = path.join(reconstruction_folder, f"{mag}.xml")
    command = f"carve --dna {input} --output {output} -i M9 -g M9 -v"
    subprocess.check_call(command, shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        metavar="INPUT",
        help="Input file (list of mags to download)",
        type=argparse.FileType("r"),
    )
    parser.add_argument("-t", "--threads", dest="n_threads", type=int, default=1)
    parser.add_argument("-o", "--output", dest="output", help="output folder")
    args = parser.parse_args()
    os.makedirs(args.output, exist_ok=True)
    with args.input as f:
        mags = f.read().splitlines()
        with ThreadPoolExecutor(max_workers=args.n_threads) as executor:
            [executor.submit(download_mag, mag, args.output) for mag in mags]
