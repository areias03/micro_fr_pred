import argparse
import micom
import pandas as pd



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        metavar="INPUT",
        help="Manifest file",
        type=str,
    )
    parser.add_argument("-o", "--output", dest="output", help="output file")
    args = parser.parse_args()
    manifest = pd.read_csv(args.input)
    with args.input as f:
