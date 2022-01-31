#!/usr/bin/env python

import gzip

# function to open both gzip'd and regular files
def gzopen(file_path, mode='rt', compresslevel=1):
    if file_path.lower().endswith('.gz'):
        return gzip.open(file_path, mode=mode, compresslevel=compresslevel)
    return open(file_path, mode)
