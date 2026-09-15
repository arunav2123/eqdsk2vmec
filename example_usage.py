#!/usr/bin/env python3
# Written by Arunav Kumar, MIT Plasma Science and Fusion center, 10th May, 2026
"""Convert a GEQDSK and optionally save the input comparison; never runs VMEC."""
import argparse
from eqdsk2vmec import convert_eqdsk_to_vmec


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('filename')
    parser.add_argument('--output', dest='output_path')
    parser.add_argument('--plot', dest='plot_path')
    parser.add_argument('--show', action='store_true')
    parser.add_argument('--mpol', type=int, default=12)
    parser.add_argument('--lasym', action='store_true', help='Preserve up-down asymmetry')
    parser.add_argument('--cocos', type=int, help='Producer COCOS: 1..8 or 11..18; omitted assumes 5 with warning')
    args = vars(parser.parse_args())
    print(convert_eqdsk_to_vmec(**args))


if __name__ == '__main__':
    main()
