'''
This script breaks down fragments of the file run.py
to carry out different tests.

A large vcf file has been create manually to test the code
'''

import os

import pandas as pd
import numpy as np

from VCF_matcher.app.run import load_sample


NAME_FILE_1 = "./test_sample.vcf"

# FIRST TEST

def test_load_sample():
    '''Verify all rows of the body of the vcf file is taken'''
    data_to_test = load_sample (NAME_FILE_1)
    assert len(data_to_test) == 10425
    return data_to_test

data_to_test = data_to_test

# SECOND TEST
def test_take_sample()
