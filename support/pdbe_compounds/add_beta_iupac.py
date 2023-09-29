"""
This script takes one PDBE CIF file and adds the `b-` prefix to the IUPAC SNFG symbol
if it is missing.
"""

import argparse
import re
import os


def setup():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=str, help="The input PDBE CIF file")
    parser.add_argument("output", type=str, help="The output PDBE CIF file")
    return parser.parse_args()


def get_conformation(line):
    """
    Return the conformation of an IUPAC CARBOHYDRATE SYMBOL line
    Either `a` or `b`
    """
    return line.split()[-1].strip()[0]


def get_snfg_symbol(line):
    """
    Return the SNFG CARBOHYDRATE SYMBOL
    """
    return line.split()[-1].strip()


def main(args):
    conf = None
    if args.input == args.output:
        args.output += ".tmp"
    with open(args.input, "r") as fi:
        with open(args.output, "w") as fo:
            for line in fi:
                if '"IUPAC CARBOHYDRATE SYMBOL"' in line:
                    conf = get_conformation(line)
                if "SNFG CARBOHYDRATE SYMBOL" in line:
                    symbol = get_snfg_symbol(line)
                    if conf == "b" and not symbol.startswith("b-"):
                        symbol = "b-" + symbol
                    line = line.replace(get_snfg_symbol(line), symbol)
                fo.write(line)
    if args.output.endswith(".tmp"):
        args.output = args.output[:-4]
        os.rename(args.output + ".tmp", args.output)
    print(f"Saved to {args.output}")


if __name__ == "__main__":
    main(setup())
    # class Args:
    #     input = "/Users/noahhk/GIT/glycosylator/support/pdbe_compounds/components.cif"
    #     output = (
    #         "/Users/noahhk/GIT/glycosylator/support/pdbe_compounds/components.beta.cif"
    #     )

    # main(Args())
