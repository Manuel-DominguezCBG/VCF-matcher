'''
This script breaks down fragments of the file run.py
to carry out different tests.

A large vcf file has been create manually to test the code
'''

import os
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
    ''' The body of the vcf file has 10425 rows'''
    assert len(data1) == 10425 


#### SECOND TEST: Does it match the sample entered by the user with the correct column?

if len(data1.columns) > 10:
    NAME_FILE_1 = "Sample_3" # Simulating that the user has entered sample 3
else:
    print('File', os.path.basename(NAME_FILE_1), 'contains one sample only.')
    NAME_FILE_1 = data1.columns.tolist()[9]

def test_2():
    '''
    The first value of the sample number 3 is "0/1:20:221,7:0.99:20:-10.2127:20"
    '''
    assert data1[NAME_FILE_1].iloc[0] == "0/1:20:221,7:0.99:20:-10.2127:20" 

#### THIRD TEST: PASS filter works ok?

# 5000 PASS variants have been put in this file (first filter)

def test_3():
    '''5000 variants I have created with FILTER == PASS'''
    assert len(data1.query('FILTER == "PASS"')) == 5000

# FOURTH TEST: VF filter works ok?

# Only 1811 variants present a VF field higher than 0.4

# To get the values from the column Sample_3
df2 = data1[NAME_FILE_1].str.split(':', expand=True)
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
    ''' 1811 variants I have created with FILTER == PASS & VF > 0.4'''
    assert len(data2) == 1811

#### FIVETH TEST: When both df (one per sample) are merge, do we miss any variant?

# To do this, I will use the same vcf file twice

# The unit identifier
data2['CHROMPOSREFALT']= data2["CHROM"].apply(str)+"."+data2["POS"].apply(str)+data2["REF"]+"."+data2["ALT"]

frames = [data2[['CHROMPOSREFALT', 'GT']],data2[['CHROMPOSREFALT', 'GT']]]

semi_final_df1 = pd.concat(frames)

def test_5():
    ''' If both dfs are merge correctly the total number of rows should be 3622'''
    assert len(semi_final_df1) == 3622

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

# Let's create two simple data frames to check that the last part of the code work as expected

data = {'CHROMPOSREFALT':['chr1.00000000001T.C', 'chr1.00000000002T.C', "chr1.00000000003T.C",'chr1.00000000004T.C'],
        'GT':["0/1", "0/0","0/0", "0/1"]}

simple_df_1 = pd.DataFrame(data, columns = ['CHROMPOSREFALT', 'GT'])
'''
        CHROMPOSREFALT   GT
0  chr1.00000000001T.C  0/1
1  chr1.00000000002T.C  0/0
2  chr1.00000000003T.C  0/0
3  chr1.00000000004T.C  0/1
'''

data = {'CHROMPOSREFALT':['chr1.00000000001T.C', 'chr1.00000000002T.C',"chr1.00000000003T.C", 'chrX.00000000003T.C'],
        'GT':["0/1", "0/0","0/1","0/1"]}

simple_df_2 = pd.DataFrame(data, columns = ['CHROMPOSREFALT', 'GT'])
'''
        CHROMPOSREFALT   GT
0  chr1.00000000001T.C  0/1
1  chr1.00000000002T.C  0/0
2  chr1.00000000003T.C  0/1
3  chrX.00000000003T.C  0/1
'''

# Bassically, these two dfs contain 1 equal heterozygous variant (row 0), 1 equal homozygous variant (row 1)
# 1 equal variant but different genotype (row 2) and 1 completely different variant.

# After this, we already know the results we should expect in the report
'''

                                                  Homozygous: 1
Number of positions with the same genotype: 2 
                                                  Heterozygous: 1
                                                
                                                  
Number of positions with different genotype: 3
                                                  


Total positions compared: 5
Percentage in common: 0.4
'''

# We can now run the last part of the code with simple_df_1 and simple_df_1
# and test if the results obtained match with the results expected.

# This is the last part of the code:

# To merge both dfs
frames = [simple_df_1[['CHROMPOSREFALT', 'GT']],simple_df_2[['CHROMPOSREFALT', 'GT']]]
semi_final_df = pd.concat(frames)

# To match variants 
final_df=(semi_final_df.assign(key=semi_final_df.groupby('CHROMPOSREFALT').cumcount())
      .pivot('CHROMPOSREFALT','key','GT')
      .rename(columns=lambda x:f"Sample{x+1}")
      .rename_axis(columns=None).reset_index())
final_df = final_df.replace(np.NaN,"99/99" ) # To avoid future error

# To compared matched variants, which are hom or het and what are complected different

final_df["Matches"] = np.where(final_df["Sample1"] == final_df["Sample2"], True, False)

    # 6.  To get hom and het 

final_df['Genotype1'] = final_df.iloc[:,2].apply(lambda x: x.split('/' or '|')[0])
final_df['Genotype2'] = final_df.iloc[:,2].apply(lambda x: x.split('/' or '|')[1])
    
final_df['Hom/het'] = np.select([final_df['Matches'] & final_df['Genotype1'].eq(final_df['Genotype2']),
                       final_df['Matches'] & final_df['Genotype1'].ne(final_df['Genotype2'])],
                       choicelist=["Hom", "Het"],
                       default=pd.NA)

# To generate the results

r4 = final_df['Matches'].value_counts().get(True, 0) # Number of positions with the same genotype
r5 = final_df['Hom/het'].value_counts().get("Hom", 0) # Number of matched homozygous variants 
r6 =final_df['Hom/het'].value_counts().get("Het", 0) # Number of matched heterozygous variants 
r7 = final_df['Matches'].value_counts().get(False, 0) # Number of positions with different genotype
r8 = len(final_df) # Total positions compared
r9 = r4/(r4+r7)    # Percentage in common

# and finally, we can run the tests

def test_6():
    '''Number of positions with the same genotype should be 2'''
    assert r4 == 2


def test_7():
    ''' Number of matched homozygous variants should 1'''
    assert r5 == 1


def test_8():
    ''' Number of matched heterozygous variants should be 1'''
    assert r6 == 1


def test_9():
    ''' Number of positions with different genotype should 3'''
    assert r7 == 3


def test_10():
    ''' Total positions compared should be 5'''
    assert r8 == 5


def test_11():
    ''' Percentage in common should be 0.4'''
    assert r9 == 0.4
