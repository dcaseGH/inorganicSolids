# quick thing to swap xyz, cif etc- may expand functionality- dhc 220217

import numpy as np


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Write input and output files mainly')
    parser.add_argument('-i', '--input', type=str, default=None, help='specify input filename ')
    parser.add_argument('-f', '--final', type=str, default=None, help='specify different filename')

    args = parser.parse_args()

    if args.input[-4:].lower() == '.xyz':
        from coreClasses import XYZFile
        structure = 
