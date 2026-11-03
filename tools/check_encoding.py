#!/usr/bin/env python

import argparse
import os

SUFFICES = [ '.f', '.F', '.f90', '.F90', '.meta' ]

parser = argparse.ArgumentParser(description='Check source files for non-ascii characters.')
parser.add_argument('--max-depth', type=int, default=None,
                    help='Maximum directory depth to search below the current directory '
                         '(0 = current directory only). Default: no limit.')
args = parser.parse_args()

top = os.getcwd()
top_depth = top.rstrip(os.sep).count(os.sep)

for root, dirs, files in os.walk(top):
    depth = root.rstrip(os.sep).count(os.sep) - top_depth
    if args.max_depth is not None and depth >= args.max_depth:
        # Don't descend any further below this directory
        dirs[:] = []
    for file in files:
        suffix = os.path.splitext(file)[1]
        if suffix in SUFFICES:
            with open(os.path.join(root, file)) as f:
                contents = f.read()
                if not contents.isascii():
                    for line in contents.split('\n'):
                        if not line.isascii():
                            raise Exception('Detected non-ascii characters in file {}, line: "{}"'.format(os.path.join(root, file), line))
