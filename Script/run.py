#!/usr/bin/env python
'''
Created on 09/12/2014

@author: manuel.dominguezbecerra@nhs.net

'''


# Import libraries 
import argparse
import pandas as pd
import os
import subprocess
import numpy as np

import dash                              # pip install dash
import dash_html_components as html
import dash_core_components as dcc
from dash.dependencies import Output, Input

from dash_extensions import Lottie       # pip install dash-extensions
import dash_bootstrap_components as dbc  # pip install dash-bootstrap-components
import plotly.express as px              # pip install plotly
from datetime import date
import calendar
from wordcloud import WordCloud          # pip install wordcloud
import sqlite3                           # pip install sqlite3
import plotly.graph_objects as go
import dash_table
from dash.dependencies import Input, Output, State


# Take the samples path/name_file

#parser = argparse.ArgumentParser()
#parser.add_argument('--vcf1', type=str, required=True)
#parser.add_argument('--vcf2', type=int, required=True)
#args = parser.parse_args()

#Name_file1 = args.vcf1 # Remove comments when finished
#Name_file2 = args.vcf2

## Delete this when finished
#Name_file1 = "/Users/monkiky/Desktop/Do-these-samples-belong-to-the-same-patient/Samples/W2008872_S8-J4PJL_copy.vcf"
Name_file1 = "/Users/monkiky/Desktop/Do-these-samples-belong-to-the-same-patient/Samples/G369884060.vcf"
Name_file2 = "/Users/monkiky/Desktop/Do-these-samples-belong-to-the-same-patient/Samples/W2110682_S9.vcf"


### 1. Remove the header of input files

# File A
cmd = "sed '/^##/ d'  {0} >  FileA.txt".format(Name_file1)
os.system(cmd)

# File B 
cmd = "sed '/^##/ d'  {0} >  FileB.txt".format(Name_file2)
os.system(cmd)

# Files to pandas object
dataA = pd.read_csv("./FileA.txt", delimiter = "\t" )
dataB = pd.read_csv("./FileB.txt", delimiter = "\t" )

# Delete files created
cmd = os. getcwd()+"/FileA.txt"
os.remove(cmd)

cmd = os. getcwd()+"/FileB.txt"
os.remove(cmd)

# If the pd dataframes have more than 10 columns
# that means there are more than 1 sample in the vcf file

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

### 2. Filter the variants  FILTER = PASS and Variant Frecuency =>0.4

# VF is in the 9th column, we need to extract and put VF data in an independent column first
df2 = dataA[Name_sample1].str.split(':', expand=True)
df2.columns = dataA.FORMAT.iloc[0].split(':')
dataA = pd.concat([dataA, df2], axis=1)


# Same for dataB
df2 = dataB[Name_sample2].str.split(':', expand=True)
df2.columns = dataB.FORMAT.iloc[0].split(':')
dataB = pd.concat([dataB, df2], axis=1)

# object to float
dataA["VF"] = dataA.VF.astype(float)
dataB["VF"] = dataB.VF.astype(float)


# ok, now we can filter and take the variant we are interested 
dataA = dataA.query('FILTER == "PASS" and VF > 0.4')
dataB = dataB.query('FILTER == "PASS" and VF > 0.4')

### 3. Let's create the Ben's Code, this is CHROM+POS+REF+ALT

dataA['CHROMPOSREDALT']=dataA["#CHROM"].apply(str)+"."+dataA["POS"].apply(str)+dataA["REF"]+"."+dataA["ALT"]

dataB['CHROMPOSREDALT']=dataB["#CHROM"].apply(str)+"."+dataB["POS"].apply(str)+dataB["REF"]+"."+dataB["ALT"]

### 4. Take common Ben's Codes in the same column and put the GT of each (both  copies) in independent columns
df3 = dataA[['CHROMPOSREDALT', 'GT']].copy()
df4 = dataB[['CHROMPOSREDALT', 'GT']].copy()
frames = [df3,df4]
semi_final_df = pd.concat(frames)



final_df=(semi_final_df.assign(key=semi_final_df.groupby('CHROMPOSREDALT').cumcount())
      .pivot('CHROMPOSREDALT','key','GT')
      .rename(columns=lambda x:f"Sample{x+1}")
      .rename_axis(columns=None).reset_index())


### 5. See how many GT match
final_df["Matches"] = np.where(final_df["Sample1"] == final_df["Sample2"], True, False)
final_df.columns = ['CHROM.POS.REF.ALT',os.path.basename(Name_file1),os.path.basename(Name_file2),"Matches" ]


### 6.1 Generate results in the terminal

# Variable to introduce in the report:
r0 = os.path.basename(Name_file1) # Get the  file name of a path/file_name 
r1 = Name_sample1
r2 = os.path.basename(Name_file2)
r3 = Name_sample2
r4 = final_df['Matches'].value_counts().get(True, 0) # Count trues if any, retunr 0
r5 = final_df['Matches'].value_counts().get(False, 0)
r6 = len(final_df)
r7 = r4/(r4+r5)
r6 == r4+r5
report = '''
 _____________________________  REPORT  ________________________________________ 

vcf 1: {0}  AND its sample name: {1}  
vcf 2:  {2}  AND its sample name: {3}

                                                  Homozigous: 
Number of positions with the same genotype: {4} 
                                                  Heterozigous:
                                                
                                                  Homozigous:
Number of positions with diferent genotype: {5} 
                                                  Heterozigous: 


Total positions compared: {6}
Percentage in common: {4}/{6}= {7}
 ____________________________ END REPORT  _______________________________________
'''
print(report.format(r0,r1,r2,r3,r4,r5,r6,r7))



print(final_df)



# draft

#(semi_final_df.pivot_table(index = ['CHROMPOSREDALT'], aggfunc ='size')==2).value_counts()[1] # position same genotipe
#position_same_genotype = (semi_final_df.pivot_table(index = ['CHROMPOSREDALT'], aggfunc ='size')==2).value_counts()[1]
#total_hom = semi_final_df[semi_final_df.GT == '1/1'].shape[0]
#total_het = semi_final_df[semi_final_df.GT != '1/1'].shape[0]