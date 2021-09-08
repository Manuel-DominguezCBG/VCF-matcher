'''
This script breaks down fragments of the file run.py 
to carry out different tests.

A large vcf file has been create manually to test the code
'''
import pytest
import os
import argparse
import pandas as pd
import numpy as np


# FIRST TEST: Check if all rows of the body of the vcf files are taken

# The script that load and save the body of the vcf into a pandas data frame object
NAME_FILE_1 = "/Users/monkiky/Desktop/VCF-matcher/Testing/test_sample.vcf"

with open(NAME_FILE_1, 'r') as f:
    for line in f:
        if line.startswith('#') and len(line)>2 and line[1] != '#':
            columns = line[1:-1].split('\t')
            break

data1 = pd.read_csv(NAME_FILE_1, comment='#', delimiter='\t', names=columns)

# The test:
def test_1():
    assert len(data1) == 10425 # The body of the vcf file has 10425 rows


#### SECOND TEST: Does it match the sample inetered by the user with the correct column?

if len(data1.columns) > 10:
    Name_sample1 = "Sample_3" # Simulating that the user has entered sample 3
else:
    print('File', os.path.basename(NAME_FILE_1), 'contains one sample only.')
    Name_sample1 = data1.columns.tolist()[9]

def test_2():
    '''
    The first value of the sample number 3 is "0/1:20:221,7:0.99:20:-10.2127:20"
    '''
    assert data1[Name_sample1].iloc[0] == "0/1:20:221,7:0.99:20:-10.2127:20" 
     
#### THIRD TEST: PASS filter works ok?

# 5000 PASS variants have been put in this file (first filter)

def test_3():
    assert len(data1.query('FILTER == "PASS"')) == 5000

# FOURTH TEST: VF filter works ok?

# Only 1811 variants present a VF fielf higher than 0.4

# To get the values from the column Sample_3
df2 = data1[Name_sample1].str.split(':', expand=True)
df2 = df2.fillna(value=np.nan)
df2 = df2.dropna(axis=1, how='any')
df2.columns = data1.FORMAT.iloc[0].split(':')
data2 = pd.concat([data1, df2], axis=1)

# Apply the filter
if "VF" in data2.columns:
    data2["VF"] = data2["VF"].astype(float)
    data2 = data2.query('VF > 0.4')

# Apply the test
def test_4():
    assert len(data2) == 1811

#### FIVETH TEST: When both df (one per sample) are merge, do we miss any variant?

# To do this, I will use the same vcf file twice

# The unit identifier
data2['CHROMPOSREFALT']= data2["CHROM"].apply(str)+"."+data2["POS"].apply(str)+data2["REF"]+"."+data2["ALT"]

frames = [data2[['CHROMPOSREFALT', 'GT']],data2[['CHROMPOSREFALT', 'GT']]]

semi_final_df = pd.concat(frames)

def test_5():
    assert len(semi_final_df) == 3622

#### SIXTH TEST: The variant matcher does what expected?

# To ensure that, it is convenient to used a two very small vcf files

# At that point of the process, a standard semi_final_df looks like this

'''
 CHROMPOSREFALT   GT
3       chr1.36933096T.C  0/1
9       chr2.25463483G.A  0/1
13      chr2.25466888G.T  1/1
16      chr2.25469502C.T  0/1
25   chr2.198267770G.GAA  0/1
..                   ...  ...
198     chrX.15838366C.T  1/1
207     chrX.39933339A.G  1/1
263    chrX.129147079T.C  1/1
264    chrX.129147373G.A  1/1
275    chrX.133511988G.A  1/1
'''

# Let's create two simplafy semi_final_df to check the part of the code that matches the variants

simple_semi_final_df_1 = 

simple_semi_final_df_2