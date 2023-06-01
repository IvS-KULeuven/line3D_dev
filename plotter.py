#!/usr/bin/env python

import line3d_utils
import os
import sys

def parse_inputarg():
    if sys.version_info.major < 3:
            raise Exception('Unsupported python version')
            quit()

    inputargs ={'Nfiles' : 0, \
                'FluxFile': []}
    if len(sys.argv[1:]) < 1:
        raise Exception('No input grid json file specified')
    else:
        argin = sys.argv.copy()
        inputargs['Nfiles'] = len(argin) - 1
        for arg in argin[1:]:
            inputargs['FluxFile'].append(arg)
        print(f"Using Flux Files {inputargs['FluxFile']}", flush = True)
    
    return inputargs


inputs = parse_inputarg()

line3d_utils.quickplot(inputs['FluxFile'])