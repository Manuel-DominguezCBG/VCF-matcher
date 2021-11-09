'''
This script breaks down fragments of the file run.py
to carry out different tests.

Input info

I have created manually a VCF file to test my application. 

This is a 24-samples VCF file with 1400 variants (rows).
In this script I open the two needed VCF files with two different approaches I have found.
the option -s of pytest allow you to introduce input so to run this test we need to type

pytest test_vcf_matcher.py -s

and then you will be ask for introducing the sample you want to used

373978487

run.py needs to VCF input, the second file is automatically take from the code
using monkeypatch (line 60). So, the same file (test_sample.vcf) is used to test application.

Then, all variants pass all filter apart from 48 variants that are LowMQ.
This is tested with test_filter_1() and test_filter_2()

Then, I test that the sample selected take the data of that column
this is done with test_take_sample()

Then, I test that the unit identifier created is what expected
this is done with test_Unit_identifier()

After this, all data is ready for comparison, so run.py created a new dataframe to combine the data needed.
to test that all PASS variants are take I test this dataframe with test_concatenate_samples()
basically ensuring that this df have two columns and the number of PASS variants X2

After this tests, I have want to check that the data transferred from samples files to the new 
comparison table match. To do this I have take 5 variants and I have ensure that the data from
sample 1 and sample 2 match with comparison table. All of this is done with test_match_variants()

From this comparison data, run.py generate the report. First thing I have done here is to ensure that
the number the hom and het match with the vcf file I created. All PASS variants are hom apart from one that is het.
This is done with test_count_hom_het_variants()

And finally I have test the rest of the values provided in the report. This is well explained in test_count_and_report()

In total, with 26 assertions I have tested all functions used in in the application.
'''

import os
import random
import pandas as pd
import numpy as np

# Import all functions from run.py
from app.run import * 

NAME_FILE_1 = "./test_sample.vcf"


def test_load_sample():
    '''Verify all rows of the body of the vcf file are taken'''
    data_to_test = load_sample (NAME_FILE_1)
    return data_to_test 
    assert len(data_to_test) == 1400
    
# To test all functions  I need the data generated in the previous functions
# That is the reason I return the output and save it.
data_to_test_1 = test_load_sample() 


def test_take_sample(monkeypatch):
    '''
    take_sample() requests the name of the column
    and take that input to split the data in new columns
    This can be tested by checking some of the first values of that
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

    #assert data_to_test_2["DP"] == "74:99"

data_to_test_2,Name_sample = take_sample(data_to_test_1,NAME_FILE_1)


def test_filter_1():
    '''
    There are 1352 variants that pass the FILTER field.
    '''
    data_to_test_3 = filter_1(data_to_test_2) # filter_1 take PASS variants
    return data_to_test_3
    assert len(data_to_test_3) == 1352

data_to_test_3 = test_filter_1()


def test_filter_2():
    '''
    The second filter takes the variants which VF is higher than 0.4
    This is a optional condition as not all vcf files present this field
    The vcf_test doesnt contains this field so the expected results would be 1352
    '''
    data_to_test_4 = filter_2(data_to_test_3) # filter_2 take VF>0.4 variants
    return data_to_test_4
    assert len(data_to_test_4) == 1352


data_to_test_4 = test_filter_2()

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
    To test this, I am going to take the vcf_test file twice and count 
    the number of columns and rows.
    Then check if the concatenation is what I expected.
    '''
    data_to_test_5 = concatenate_samples(data_to_test_4,data_to_test_4)
    return data_to_test_5
    assert data_to_test_5.shape[0] == 2       # Two columns
    assert data_to_test_5.shape[1] == 1352*2  # Double of rows

data_to_test_5 = test_concatenate_samples()


def test_match_variants():
    '''
    I have picked up a few  variants from the original files
    and I have done this variants manuallyI can check now if there these matches 
    1.12026355A.C   0/1
    3.52410008G.C   0/0
    7.810263C.T     0/0
    22.26860269G.C  1/1
    X.38145690TCCTTCCTCCTCTTCCCCCTCC.T 0/0
    '''
    data_to_test_6 = match_variants(data_to_test_5)

    return data_to_test_6

    assert data_to_test_6.loc[data_to_test_6['CHROMPOSREFALT'] == "1.12026355A.C", 'Sample1'].item() == "0/1"
    assert data_to_test_6.loc[data_to_test_6['CHROMPOSREFALT'] == "1.12026355A.C", 'Sample2'].item() == "0/1"

    assert data_to_test_6.loc[data_to_test_6['CHROMPOSREFALT'] == "3.52410008G.C", 'Sample1'].item() == "0/0"
    assert data_to_test_6.loc[data_to_test_6['CHROMPOSREFALT'] == "3.52410008G.C", 'Sample2'].item() == "0/0"

    assert data_to_test_6.loc[data_to_test_6['CHROMPOSREFALT'] == "7.810263C.T", 'Sample1'].item() == "0/0"
    assert data_to_test_6.loc[data_to_test_6['CHROMPOSREFALT'] == "7.810263C.T", 'Sample2'].item() == "0/0"

    assert data_to_test_6.loc[data_to_test_6['CHROMPOSREFALT'] == "22.26860269G.C", 'Sample1'].item() == "1/1"
    assert data_to_test_6.loc[data_to_test_6['CHROMPOSREFALT'] == "22.26860269G.C", 'Sample2'].item() == "1/1"

    assert data_to_test_6.loc[data_to_test_6['CHROMPOSREFALT'] == "X.38145690TCCTTCCTCCTCTTCCCCCTCC.T", 'Sample1'].item() == "0/0"
    assert data_to_test_6.loc[data_to_test_6['CHROMPOSREFALT'] == "X.38145690TCCTTCCTCCTCTTCCCCCTCC.T", 'Sample2'].item() == "0/0"

data_to_test_7 = test_match_variants()
data_to_test_7["Matches"] = np.where(data_to_test_7["Sample1"] == data_to_test_7["Sample2"], True, False)
data_to_test_7.columns = ['CHROMPOSREFALT',os.path.basename(NAME_FILE_1),os.path.basename(NAME_FILE_1),"Matches" ]

def test_count_hom_het_variants():
    '''
    In sample 373978487 there is only one heterozygous, the rest of the variants are homozygous
    I have done this manually
    ''' 
    data_to_test_7['Hom/het'] = count_hom_het_variants(data_to_test_7)

    assert data_to_test_7['Hom/het'].value_counts().get("Hom", 0) == 1341
    assert data_to_test_7['Hom/het'].value_counts().get("Het", 0) == 1
    return data_to_test_7


data_to_test_8 = test_count_hom_het_variants()

def test_count_and_report():
    
    r4,r5,r6,r7,r8,r9 = count_and_report(data_to_test_7)
    assert r4 == 1342 # because we are testing the same vcf file twice
                      # the number of variant with the same position is equal to the total
                      # number of variants 
    assert r5 == 1341 # 1342 - 1 heterozygous
    assert r6 == 1    # I have create only one heterozygous
    assert r7 == 0    # 0 positions with different genotype because the samples compared are the same
    assert r8 == 1342 # Total number of variants
    assert r9 == 1.0  # Percentage in common is 100% because the samples compared are the same.




# Here I test the VF filter using a second vcf file I have created.



NAME_FILE_2 = "./to_test_VF.vcf"

data_VF = load_sample(NAME_FILE_2)
data_VF,Name_sample_1 = take_sample(data_VF,NAME_FILE_2)


def test_filter_3():
    '''
    There are 8 variants in to_test_VF.vcf.
    6 PASS VF>0.4
    1 PASS VF<0.4
    1 lowMQ VF>0.4
    So, the filter_1 should report 7 because there is only one lowMQ variant
    and filter_2 should report 6  (8-1-1 =6)
    '''
    PASS_variants = filter_1(data_VF)
    assert len(PASS_variants) == 7

    PASS_FV_variants = filter_2(PASS_variants)
    assert len(PASS_FV_variants) == 6

