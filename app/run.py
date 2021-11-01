#!/usr/bin/env python
'''
Created on 07/09/2021
@author: manuel.dominguezbecerra@nhs.net
'''

# Import libraries
import os
import argparse
import pandas as pd
import numpy as np
import sys

### 1. Load the files

# Take only the data contains in the body of the vcf files.

def load_sample (Name_file):
    '''
    Take the header and the body of the CSV file. Ignore metada.
    '''
    with open(Name_file, 'r') as f:
        for line in f:
            if line.startswith('#') and len(line)>2 and line[1] != '#':
                columns = line[1:-1].split('\t')
                data = pd.read_csv(Name_file, comment='#', delimiter='\t', names=columns)
                break
    return data

    ### 2. One or > one sample in the files? 

def take_sample(data,NAME_FILE):
    '''
    VCF files can have 1 or > 1 samples
    We need to know what sample the user wants to compare
    After this, the data of the sample is all stored in one  column
    to access the data, it needs to  split it in different columns
    '''
    if len(data.columns) > 10:
        Name_sample = input('It seems that {} has more than one sample. Please introduce the name of the sample you wish to analyse:  '.format(os.path.basename(NAME_FILE)))
    else:
        print('File', os.path.basename(NAME_FILE), 'contains one sample only.')
        Name_sample = data.columns.tolist()[9]
    
    tmp_df_2 = data[Name_sample].str.split(':', expand=True)
    # I have seen very unlikely Na values
    # To solve this, they are deleted
    tmp_df_2 = tmp_df_2.fillna(value=np.nan)
    tmp_df_2 = tmp_df_2.dropna(axis=1, how='any')
    # To get the name of the new columns created
    tmp_df_2.columns = data.FORMAT.iloc[0].split(':')
    return pd.concat([data, tmp_df_2], axis=1),Name_sample
    
    ### 3. Filter the variants

def filter_1(data):
    '''
    Take only the variants that "PASS" the "FILTER" field
    '''
    return data.query('FILTER == "PASS"')


def filter_2(data):
    '''
    Take the variants which VF > 0.4
    '''
    if "VF" in data.columns:
        data["VF"] = data["VF"].astype(float)
        return data.query('VF > 0.4')
    else:
        return data

   ### 4. Variants comparison

def Unit_identifier(data):
    '''
    Two variants are the same if  CHROM, POS, REF, ALT and GT are equal
    Here, the unit identifier is created
    '''
    return data["CHROM"].apply(str)+"."+data["POS"].apply(str)+data["REF"]+"."+data["ALT"]


def concatenate_samples(data1,data2):
    '''
    To proceeded, both data frames need to be concatenate 
    '''
    frames = [data1[['CHROMPOSREFALT', 'GT']],data2[['CHROMPOSREFALT', 'GT']]]
    return pd.concat(frames)

def match_variants(data):
    '''
    Match common variants in one column and put the GT values of both samples
    in different columns.
    No common variants are also set in the first column giving NaN values 
    in the third column. I have solved this by substituting the NaN for 99/99 values
    '''
    matched_df = (data.assign(key=data.groupby('CHROMPOSREFALT').cumcount())
        .pivot('CHROMPOSREFALT','key','GT')
        .rename(columns=lambda x:f"Sample{x+1}")
        .rename_axis(columns=None).reset_index())
    matched_df = matched_df.replace(np.NaN,"99/99" )
    return matched_df.iloc[:, : 3]

def count_hom_het_variants(data):
    data['Genotype1'] = data.iloc[:,2].apply(lambda x: x.split('/' or '|')[0])
    data['Genotype2'] = data.iloc[:,2].apply(lambda x: x.split('/' or '|')[1])
    return np.select([data['Matches'] & data['Genotype1'].eq(data['Genotype2']),
                       data['Matches'] & data['Genotype1'].ne(data['Genotype2'])],
                       choicelist=["Hom", "Het"],
                       default=pd.NA)

def count_and_report(data):
    r4 = data['Matches'].value_counts().get(True, 0) # Count trues if any, return 0
    r5 = data['Hom/het'].value_counts().get("Hom", 0)
    r6 = data['Hom/het'].value_counts().get("Het", 0)
    r7 = data['Matches'].value_counts().get(False, 0)
    r8 = len(data)
    r9 = r4/(r4+r7)
    return r4,r5,r6,r7,r8,r9 


if __name__ == "__main__":
    # execute only if run as a script

    # To allow the user to introduce the file after running "python run.ppy"

    # Example:    python run.py --vcf1 /path/file_name1 --vcf2 /path/file_name2

    # Enter the path/file names

    # Create the parser
    parser = argparse.ArgumentParser(description="Compare 2 VCF files and return the number of common variants that genotype match. Full documentation here: https://github.com/Manuel-DominguezCBG/VCF-matcher")

    # Add the arguments
    parser.add_argument('vcf1',metavar="vcf1", type=str, help="First VCF input. Eg.: PATH/filename.vcf")
    parser.add_argument('vcf2',metavar="vcf2", type=str, help="Second VCF input, Eg.: PATH/filename.vcf")
 
    # Execute the parse_args() methods
    args = parser.parse_args()
    NAME_FILE_1 = args.vcf1 
    NAME_FILE_2 = args.vcf2


    # Load the files
    dataA = load_sample(NAME_FILE_1)
    dataB = load_sample(NAME_FILE_2)

    # Save now the data of this sample and its name (we will need for the report)
    dataA,Name_sample_1 = take_sample(dataA,NAME_FILE_1)
    dataB,Name_sample_2 = take_sample(dataB,NAME_FILE_2)

    # Filter 1
    dataA = filter_1(dataA)
    dataB = filter_1(dataB)

    # Filter 2
    dataA = filter_2(dataA)
    dataB = filter_2(dataB)

    # Unit_identifier
    dataA['CHROMPOSREFALT'] = Unit_identifier(dataA)
    dataB['CHROMPOSREFALT'] = Unit_identifier(dataB)

    # Concatenate both data frames
    semi_final_df = concatenate_samples(dataA,dataB)

    # Match common variants
    final_df= match_variants(semi_final_df)

#### Generate report ####

    # Variables needed for both type of reports

    r0 = os.path.basename(NAME_FILE_1) # Get the  file name of a path/file_name
    r1 = Name_sample_1                 # Get the name of the sample 1
    r2 = os.path.basename(NAME_FILE_2)
    r3 = Name_sample_2

    if "Sample2" not in final_df.columns:
    # If Sample2 is not in final_df, that means there is not common position between the samples.
    # This is a the most simple kind of report with only the name of the files
    # and samples and informing that no common positions have been found between the samples
        REPORT_A = '''
 _____________________________  REPORT  ________________________________________ 

   vcf 1: {0}  AND its sample name: {1}  
   vcf 2:  {2}  AND its sample name: {3}

   No common positions found between samples

 ____________________________ END REPORT  _______________________________________
'''
        print(REPORT_A.format(r0,r1,r2,r3,))

    # If Sample2 is in final_df.columns, that means that at least one common variant has been found
    # so, generate the data from these common variants
    else:
    # Else count common positions and carry on the report.
        final_df["Matches"] = np.where(final_df["Sample1"] == final_df["Sample2"], True, False)
        final_df.columns = ['CHROMPOSREFALT',os.path.basename(NAME_FILE_1),os.path.basename(NAME_FILE_2),"Matches" ]
        print(final_df)
    #  Get and count hom and het common variants
        final_df['Hom/het'] = count_hom_het_variants(final_df)

    # Generate the second type of report
        r4,r5,r6,r7,r8,r9 = count_and_report(final_df)
        REPORT_B = '''
 _____________________________  REPORT  ________________________________________ 

vcf 1: {0}  AND its sample name: {1}  
vcf 2:  {2}  AND its sample name: {3}

                                                  Homozygous: {5} 
Number of positions with the same genotype: {4} 
                                                  Heterozygous: {6} 
                                                
                                                  
Number of positions with different genotype: {7} 
                                                  


Total positions compared: {8}
Percentage in common: {4}/{8}= {9}

 ____________________________ END REPORT  _______________________________________
'''
        print(REPORT_B.format(r0,r1,r2,r3,r4,r5,r6,r7,r8,r9))

