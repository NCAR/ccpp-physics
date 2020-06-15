#!/usr/bin/env python

#import chardet
import os
import sys


SUFFICES = [ '.f', '.F', '.f90', '.F90', '.meta' ]

for root, dirs, files in os.walk(os.getcwd()):
    #print root, dirs, files
    for file in files:
        suffix = os.path.splitext(file)[1]
        #print file, suffix
        if suffix in SUFFICES:
            with open(os.path.join(root, file)) as f:
                contents = f.read()
            try:
                contents.decode('ascii')
            except UnicodeDecodeError:
                for line in contents.split('\n'):
                    try:
                        line.decode('ascii')
                    except UnicodeDecodeError:
                        raise Exception('Detected non-ascii characters in file {}, line: "{}"'.format(os.path.join(root, file), line))
