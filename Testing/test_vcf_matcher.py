'''
This script breaks down fragments of the file run.py
to carry out different tests.

A large vcf file has been create manually to test the code
'''

import os
import random
import pandas as pd
import numpy as np

# Import all functions from run.py
from VCF_matcher.app.run import * 


NAME_FILE_1 = "./test_sample.vcf"

# FIRST TEST
def test_load_sample():
    '''Verify all rows of the body of the vcf file are taken'''
    data_to_test = load_sample (NAME_FILE_1)
    return data_to_test # To concadenate the test I need the data of the functions
    assert len(data_to_test) == 1400
    
# To concadenate the test I need the data of the functions
data_to_test_1 = test_load_sample() 

# SECOND TEST
def test_take_sample(monkeypatch):
    '''
    take_sample() requests the name of the column
    and take that input to split the data in new columns
    This can be tested by checking some ofthe first values of that
    columns (GT:AD:DP:GQ:PL)
    monkeypatch simulate the input of the user
    Another option would be by introducing the input manually by writing 
    pytest test_vcf_matcher.py -s
    '''
    monkeypatch.setattr('builtins.input', lambda _: "373978487") # The 9th sample in the file
    data_to_test_2 = take_sample(data_to_test_1,NAME_FILE_1)
    return data_to_test_2
    assert data_to_test_2["GT"] == "0/1"
    assert data_to_test_2["AD"] == "28,46"
    assert data_to_test_2["DP"] == "74:99"

data_to_test_2,Name_sample = take_sample(data_to_test_1,NAME_FILE_1)

# THIRD TEST
def test_filter_1():
    '''
    There are 1352 variants that pass the FILTER field.
    '''
    data_to_test_3 = filter_1(data_to_test_2)
    return data_to_test_3
    assert len(data_to_test_3) == 1352

data_to_test_3 = test_filter_1()
# 4TH TEST
def test_filter_2():
    '''
    The second filter takes the variants which VF is higher than 0.4
    This is a optional condition as not all vcf files present this field
    The vcf_test doesnt contains this field so the expected results would be 1352
    '''
    data_to_test_4 = filter_2(data_to_test_3)
    return data_to_test_4
    assert len(data_to_test_4) == 1352


data_to_test_4 = test_filter_2()
# 5TH TEST

def test_Unit_identifier():
    '''
    Unit_identifier() creates a new columns concatenating CHROM,POS,REF and ALT
    We can test this matching that the concatenation of these columns is equal that 
    the new columns generate by this function in one ramdon row
    '''
    num1 = random.randint(0, len(data_to_test_4))
    data_to_test_5 = Unit_identifier(data_to_test_4)
    return data_to_test_5
    CHROM = data_to_test_5["CHROM"].iloc[num1]
    POS = data_to_test_5["POS"].iloc[num1]
    REF = data_to_test_5["REF"].iloc[num1]
    ALT = data_to_test_5["ALT"].iloc[num1]
    CHROMPOSREFALTT = data_to_test_5["CHROMPOSREFALTT"].iloc[num1]
    assert CHROMPOSREFALTT == CHROM + "." + POS + REF + "." + ALT

data_to_test_4['CHROMPOSREFALT'] = test_Unit_identifier()


def test_concatenate_samples():
    '''
    concatenate_samples concatenates the columns CHROMPOSREFALT and GT of both files.
    To test this, I am goint to take the vcf_test file twice and count 
    the number of columns and rows.
    Then check if the concatenation is what I expected.
    '''
    data_to_test_5 = concatenate_samples(data_to_test_4,data_to_test_4)
    return data_to_test_5
    assert data_to_test_5.shape[0] == 2       # Two columns
    assert data_to_test_5.shape[1] == 1352*2  # Double of rows

data_to_test_5 = test_concatenate_samples()

# 373978487

print(match_variants(data_to_test_5))

def test_match_variants():
    '''
    I have picked up a couple of variants from the original file
    and I have done this variants manuallyI can check now if there these matches 
    1.12026355A.C   0/1
    3.52410008G.C   0/0
    7.810263C.T     0/0
    22.26860269G.C  0/1
    X.38145690TCCTTCCTCCTCTTCCCCCTCC.T 0/0
    '''
    data_to_test_6 = match_variants(data_to_test_5)
    assert data_to_test_6.loc[data_to_test_6['CHROMPOSREFALT'] == "1.12026355A.C", 'Sample1'].item() == "0/1"
    assert data_to_test_6.loc[data_to_test_6['CHROMPOSREFALT'] == "1.12026355A.C", 'Sample2'].item() == "0/1"