from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from actsims import noise,utils
from soapack import interfaces as sints
import numpy as np
import os,sys
from enlib import bench
import argparse

# Parse command line
parser = argparse.ArgumentParser(description='Make covsqrt, generate some test sims, make verification plots.')
parser.add_argument("version", type=str,help='A prefix for a unique version name')
parser.add_argument("model", type=str,help='Name of a datamodel specified in soapack.interfaces.')
parser.add_argument("--mask-kind", type=str,  default="binary_apod",help='Mask kind')
parser.add_argument("--mask-patch", type=str,  default=None,help='Mask patch')
parser.add_argument("--mask-pad", type=int,  default=None,
                    help='Mask additional padding. No padding is applied to the extracted mask if any specified.')
parser.add_argument("--extract-mask", type=str,  default=None,
                    help='Make sims on the big mask but do all the analysis on an extract of this version.')
parser.add_argument("--binary-percentile", type=float,  default=10.,help='Binary percentile for sim masking.')
parser.add_argument("--season", type=str,help='Season')
parser.add_argument("--array", type=str,help='Array')
parser.add_argument("--patch", type=str,help='Patch')
args = parser.parse_args()


mask_patch = args.patch if args.mask_patch is None else args.mask_patch

with bench.show("region"):
    if args.extract_mask is not None:
        region_map = sints.get_act_mr3_crosslinked_mask(mask_patch,
                                                        version=args.extract_mask,
                                                        kind=args.mask_kind,
                                                        season=args.season)
    else:
        region_map = None

with bench.show("ngen init"):
    ngen = noise.NoiseGen(args.version,extract_region=region_map,ncache=1)

with bench.show("make a sim"):
    ngen.generate_sim(season=args.season,patch=args.patch,array=args.array,seed=1)

with bench.show("make another sim"):
    ngen.generate_sim(season=args.season,patch=args.patch,array=args.array,seed=2)
