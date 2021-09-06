# VCF-matcher

Determine if two VCF files came from the same biological source or patient by comparing common variants.

## Introduction

In facilities that manage high sample throughput (E.g. in a clinical genomics lab) sample mix-up or mislabelling is a common problem ( Koboldt et al. , 2010; Grimm et al. , 2010 ). This can lead to incorrect data processing and analysis which conflicting results and error conclusions/reports might have a huge negative impact in the diagnostic provided for the clinical genomics lab in which I am working now.

In one of the discussions carried out with my supervisor, he showed me an approach of determining manually whether two Variant Call Format (VCF) files represent samples from the same biological source by comparing their genotypes. Following this method, I propose here a script of Python that performs exactly the same, it provides a rapid pair-wise comparison of two VCF files.

##  Implementation, features and installation

This program is a Python command-line tool (Python v3.9) for Linux operating systems. The script only requires two VCF files and 4 Python libraries. The Python libraries are:

- argparse
- pandas
- os
- numpy

For installation, I suggest cloning this repository to your machine and set up a virtual environment so that it doesn't mess up any other python installation you've got. This is very simple to do by following the instructions [here](https://docs.python.org/3/tutorial/venv.html). A very condensed version:

0. In a folder, e.g.  `/home/Name/NGS` 

1. Create the environment with a command like this:  
    `python3 -m venv VCF-matcher`
Or using Conda:
    `conda create --prefix VCF-matcher`

2. Activate the environment.
    `cd /home/Name/NGS/VCF-matcher/bin/activate`
Or using Conda:
    `conda activate VCF-matcher`

This activates the virtual environment. You will need to do it every time you want to work on your code.

3. Install the required python modules. To do this, you can use the requirements.txt file from the repository. Something like:
    ` pip install -r requirements.txt`

4. Run the test script to verify that everything is installed and ready to go:
    `python test_installation.py`

5. Run the script
    `python run.py --vcf1 Path/file_name_of_sample1 --vcf2 Path/file_name_of_sample2`


## Functional/non-functional requirements

Detailed explanations is provided in the code. In the next 2 UML diagrams is in general the process.

Diagram 1

1. 2 VCF files as input
2. Programm takes the body of the files (ignores meta-data lines) and check the number columns. The body of a VCF (Version 4.2) is formed by the next mandatory columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO and FORMAT. Then, it is followed by the data of the samples. More than 10 columns means that there are more than 1 sample so the program ask the user to select the desirable sample. 10 columns means that it is a one-sample VCF.
3. Creates 2 data frames (one per sample) with the data we need.
4. Variant filtering. A couple of filters to ignore somactics variants or variants with low quality.

Diagram 2 

5. Generate a unit identifier. This is the concatenation of chromosome, position, reference and alteration alleles. It also takes the genotypes (GT) values.
6. Merge both data frames and count common variants (variants that have the same unit identifier). Count the number of variants with equal unit identifier and GT. Finally, count the number of heterozigous and homozigous common variants.
7. Generate the report. The different type of reports are explain below.


Every report starts with two lines informing the name of the VCF files and the samples.

After this, the number of positions with the same genotype ( unit of identifier and GT) and the number of homozygous and heterozygous. Then, the sample of positions with a different genotype. Finally, the proportion of variants that match. 

When both vcf files are the same.
( File names and sample names in these results have been modified to avoid any data protection issue.)
```
 _____________________________  REPORT  ________________________________________ 

vcf 1: 1  AND its sample name: A 
vcf 2:  1  AND its sample name: A
                                                  Homozigous: 20 
Number of positions with the same genotype: 62 
                                                  Heterozigous: 42 
                                                
Number of positions with different genotype: 0 
                                                  
Total positions compared: 62
Percentage in common: 62/62= 1.0
 ____________________________ END REPORT  _______________________________________
```

Typical results when both samples belong to the same patient

```

 _____________________________  REPORT  ________________________________________ 

vcf 1: 2.vcf  AND its sample name: AAA  
vcf 2:  3.vcf  AND its sample name: BBB
                                                  Homozigous: 20 
Number of positions with the same genotype: 59 
                                                  Heterozigous: 39 
Number of positions with different genotype: 4 
                                                  
Total positions compared: 63
Percentage in common: 59/63= 0.9365079365079365
 ____________________________ END REPORT  _______________________________________
```

Typical results when  samples don´t belong to the same patient

```
File 12.vcf contains one sample only.
File 20.vcf contains one sample only.

 _____________________________  REPORT  ________________________________________ 

vcf 1: 12.vcf  AND its sample name: SSS  
vcf 2:  20.vcf  AND its sample name: XXX
                                                  Homozigous: 19 
Number of positions with the same genotype: 33 
                                                  Heterozigous: 14 
                                                
Number of positions with different genotype: 45 
                                                  
Total positions compared: 78
Percentage in common: 33/78= 0.4230769230769231
 ____________________________ END REPORT  _______________________________________
```

When there is not any common variant the report  looks like this 

```
 _____________________________  REPORT  ________________________________________ 

vcf 1:   AND its sample name: 
vcf 2:    AND its sample name: 

   No common positions found between samples

 ____________________________ END REPORT  _______________________________________
```
This last case is an indicative of the vcf files represent different kind of test (eg. Myeloid_1.2 and genotyping). In this cirscuntances, even samples from the same patient will show 0 positions in common. 



##  A small experiment to check if the methodology follows works
In a small experiment carried out with 40 samples from the same group of patients versus  X samples from different patients. I have been able to show how the fraction of sites with a common genotype ranges from 80 to 100%. In samples belonging to different patients, this value falls below 50%. 


After running the script with samples belonging to the same patient (matched samples) and samples from different patients (unmatched samples) I have compared the results in the following scatter plot. This shows the results of X pair-wise comparisons between 40 samples belonging to the same set of individuals (matched samples)  and X samples not belonging to the same individuals (unmatched samples). It can be seen how samples that belong to the same biological source present a higher proportion of positions with common genotypes (X-axis).

(Scatter plot here)