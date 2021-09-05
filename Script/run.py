#!/usr/bin/env python
'''

Created on 

@author: manuel.dominguezbecerra@nhs.net

'''


# Import libraries 
import argparse                          # pip install argparse
import pandas as pd                      # pip install dash
import os
import numpy as np




# To allow the user to introduce the file after running "python run.ppy"

# Example:    python run.py -- /path/file_name1 -- /path/file_name2

   ### 1. Take the files
    
# During the development of this application, I have needed to run the script many times
# To do this in a easy way, I directly introduce here the path/file_name in a variable
# As it can be seen below

Name_file1 = "/Users/monkiky/Desktop/VCF-matcher/Samples/Same_patients/1/W2103070_S8.vcf"

Name_file2 = "/Users/monkiky/Desktop/VCF-matcher/Samples/Same_patients/1/W1815770_S16.vcf"

# In the script, the previous two lines doesn´t appear and the following lines are incommeneted

#parser = argparse.ArgumentParser()
#parser.add_argument('--vcf1', type=str, required=True)
#parser.add_argument('--vcf2', type=int, required=True)
#args = parser.parse_args()
#Name_file1 = args.vcf1 
#Name_file2 = args.vcf2

# This 6 lines allow the user to select the two neccesary files to run the script

#During the development, the few lines make running the script easier.
# In the run.py script, the previous lines are uncommented and the following lines are deleted.




# From here to the end the following lines are identical that the one find in the script


   ### 2. Take only the data contains in the vcf files.

# Originally what I did was to delete the metedata and save in a new file the data
# This approach was over complicated and doesn´t work in window machines.
# Instead of this, the new approach directly take the data and save it in pandas dataframe object.
# This is also a more efficient method.

# For the file number 1
with open(Name_file1, 'r') as f:
    for line in f:
        if line.startswith('#') and len(line)>2 and line[1] != '#':
            columns = line[1:-1].split('\t')
            break

dataA = pd.read_csv(Name_file1, comment='#', delimiter='\t', names=columns)

# For the second file
with open(Name_file2, 'r') as f:
    for line in f:
        if line.startswith('#') and len(line)>2 and line[1] != '#':
            columns = line[1:-1].split('\t')
            break
            
dataB = pd.read_csv(Name_file2, comment='#', delimiter='\t', names=columns)


# If the pd data frame has more than 10 columns
# that means there are more than 1 sample in the vcf file
# If more than one file is found, next lines ask which sample to take.

if len(dataA.columns) > 10:
     Name_sample1 = input('It seems that {} has more than one sample. Please introduce the name of the sample you wish to analyse:  '.format(os.path.basename(Name_file1)))
else:
    print('File', os.path.basename(Name_file1), 'contains one sample only.')
    Name_sample1 = dataA.columns.tolist()[9]

if len(dataB.columns) > 10:
     Name_sample2 = input('It seems that {} has more than one sample. Please introduce the name of the sample you wish to analyse:  '.format(os.path.basename(Name_file2)))
else:
    print('File', os.path.basename(Name_file2), 'contains one sample only.')
    Name_sample2 = dataB.columns.tolist()[9]

    ### 3. Filter the variants  
    
    
# Two filters are applied. First, I take only the variants that "PASS" the "FILTER" field. 
# This is a mandatory field.

# Then, I take only the variants which Variant Frecuency (VF) is highter or eqqual than. 0.4
# This is not a mandatory field in vcf filed

# VF is in the 9th column with another data. Example: GT:GQ:AD:VF:NL:SB:GQX   0/1:100:3501,254:0.0676:20:-100.0000:100
#we need to extract and put VF data in an independent column first
df2 = dataA[Name_sample1].str.split(':', expand=True)

df2 = df2.fillna(value=np.nan)
df2 = df2.dropna(axis=1, how='any')

df2.columns = dataA.FORMAT.iloc[0].split(':')
dataA = pd.concat([dataA, df2], axis=1)

# Same for dataB
df3 = dataB[Name_sample2].str.split(':', expand=True)
df3 = df3.fillna(value=np.nan) 
df3 = df3.dropna(axis=1, how='any') 

df3.columns = dataB.FORMAT.iloc[0].split(':')
dataB = pd.concat([dataB, df3], axis=1)

# Now, let´s filter the variant
# We dont want variants that doesnt pass the Filter status
# This field is mandatory, no condition is needed because it is applied to any type of VCF file.
dataA = dataA.query('FILTER == "PASS"')
dataB = dataB.query('FILTER == "PASS"')
                    
# Now, I apply Variant Frecuency (VF) filter (non-mandatory field)

# If there is a VF Genotype fields in the file
# Take the variants that are > 0.4
# If not, no VF filter is applied

if "VF" in dataA.columns:
    dataA["VF"] = dataA["VF"].astype(float)
    dataA = dataA.query('VF > 0.4')
    
if "VF" in dataB.columns:
    dataB["VF"] = dataB["VF"].astype(float)
    dataB = dataB.query('VF > 0.4')
    

   ### 4. Variants comparison

# Two variants are the same if  CHROM, POS, REF, ALT and GT are equal

dataA['CHROMPOSREDALT']= dataA["CHROM"].apply(str)+"."+dataA["POS"].apply(str)+dataA["REF"]+"."+dataA["ALT"]

dataB['CHROMPOSREDALT']= dataB["CHROM"].apply(str)+"."+dataB["POS"].apply(str)+dataB["REF"]+"."+dataB["ALT"]


df3 = dataA[['CHROMPOSREDALT', 'GT']].copy()
df4 = dataB[['CHROMPOSREDALT', 'GT']].copy()
frames = [df3,df4]
semi_final_df = pd.concat(frames)
# semi_final_df saves the CHROM, POS, REF, ALT in the first column and the GT in the second column


#Example:
'''
CHROMPOSREDALT     GT

chr1.123456T.C     0/1       From sample 1
chr2.123456A.C     1/1       From sample 1
chr1.123456T.C     0/1       From sample 2

'''


# Let´s reduce the number of rows
final_df=(semi_final_df.assign(key=semi_final_df.groupby('CHROMPOSREDALT').cumcount())
      .pivot('CHROMPOSREDALT','key','GT')
      .rename(columns=lambda x:f"Sample{x+1}")
      .rename_axis(columns=None).reset_index())

# Final_df contains all common and non-common variants that passed the filters

# The pivot command gives NaNs in the column of the second sample when a variant is only found in one of the samples
# This NaN values return errors in the following lines. This is an error I have found when the script was developed
# To solve the issue without introducing many changes, NaN values are replace for 99/99 values.
final_df = final_df.replace(np.NaN,"99/99" )

   ### 4. See how many GT match

# Variables needed for both type of reports

r0 = os.path.basename(Name_file1) # Get the  file name of a path/file_name 
r1 = Name_sample1                 # Get the name of the sample 1
r2 = os.path.basename(Name_file2)
r3 = Name_sample2

if "Sample2" not in final_df.columns:
    # If Sample2 is not in final_df, that means there is not common position between the samples.
    report0 = '''
 _____________________________  REPORT  ________________________________________ 

vcf 1: {0}  AND its sample name: {1}  
vcf 2:  {2}  AND its sample name: {3}

   No common positions found between samples

 ____________________________ END REPORT  _______________________________________
'''
    print(report0.format(r0,r1,r2,r3,))

else:    
    # else match common position and carry on the report.
    final_df["Matches"] = np.where(final_df["Sample1"] == final_df["Sample2"], True, False)
    final_df.columns = ['CHROM.POS.REF.ALT',os.path.basename(Name_file1),os.path.basename(Name_file2),"Matches" ]

    # 6.  To get hom and het 

    final_df['Genotype1'] = final_df.iloc[:,2].apply(lambda x: x.split('/' or '|')[0])
    final_df['Genotype2'] = final_df.iloc[:,2].apply(lambda x: x.split('/' or '|')[1])
    
    final_df['final_match'] = np.select([final_df['Matches'] & final_df['Genotype1'].eq(final_df['Genotype2']),
                       final_df['Matches'] & final_df['Genotype1'].ne(final_df['Genotype2'])],
                       choicelist=[True, False],
                       default=pd.NA)

   ### 6. Generate results in the terminal

# Variable to introduce in the report :
    
    r4 = final_df['Matches'].value_counts().get(True, 0) # Count trues if any, retunr 0
    r5 = final_df['final_match'].value_counts().get(True, 0)
    r6 =final_df['final_match'].value_counts().get(False, 0)
    r7 = final_df['Matches'].value_counts().get(False, 0)
    r8 = len(final_df)
    r9 = r4/(r4+r7)
    report = '''
 _____________________________  REPORT  ________________________________________ 

vcf 1: {0}  AND its sample name: {1}  
vcf 2:  {2}  AND its sample name: {3}

                                                  Homozigous: {5} 
Number of positions with the same genotype: {4} 
                                                  Heterozigous: {6} 
                                                
                                                  
Number of positions with different genotype: {7} 
                                                  


Total positions compared: {8}
Percentage in common: {4}/{8}= {9}

 ____________________________ END REPORT  _______________________________________
'''
    print(report.format(r0,r1,r2,r3,r4,r5,r6,r7,r8,r9))
    
