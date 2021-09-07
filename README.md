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
    `conda create --name vcf_matcher  python --no-default-packages`

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


## Functional/non-functional requirements and results

Detailed explanations is provided in the code. However, the general process is explained here.

### Diagram 1

![alt text](https://github.com/Manuel-DominguezCBG/VCF-matcher/blob/main/Images_slides_and_stuff_to_explain_the_application/Slide1.JPG?raw=true)

### Diagram 2

![alt text](https://github.com/Manuel-DominguezCBG/VCF-matcher/blob/main/Images_slides_and_stuff_to_explain_the_application/Slide2.JPG?raw=true)


Diagram 1 explanation

1. 2 VCF files as input
2. Program takes the body of the files (ignores meta-data lines) and check the number columns. The body of a VCF (Version 4.2) is formed by the next mandatory columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO and FORMAT. Then, it is followed by the data of the samples. More than 10 columns means that there are more than 1 sample so the program ask the user to select the desirable sample. 10 columns means that it is a one-sample VCF.
3. Creates 2 data frames (one per sample) with the data we need.
4. Variant filtering. A couple of filters to ignore somactics variants or variants with low quality.

Diagram 2 explanation

5. Generate a unit identifier. This is the concatenation of chromosome, position, reference and alteration alleles. It also takes the genotypes (GT) values.
6. Merge both data frames and count common variants (variants that have the same unit identifier). Count the number of variants with equal unit identifier and GT. Finally, count the number of heterozigous and homozigous common variants.
7. Generate the report. The different type of reports are explain below.


Every report starts always with two lines informing the name of the VCF files and the samples.

After this, the number of positions with the same genotype ( unit of identifier and GT) and the number of homozygous and heterozygous. Then, the sample of positions with a different genotype. Finally, the proportion of variants that match. I have copied and page some reports under different and typical cirscunstances.

#### When both vcf files are the same.

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

#### When both samples belong to the same patient

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

#### When  samples don´t belong to the same patient

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

#### When there is not any common variant the report looks like this 

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

![alt text](https://github.com/Manuel-DominguezCBG/VCF-matcher/blob/main/Images_slides_and_stuff_to_explain_the_application/download.png?raw=true)




The following documentation is not necessary to undertand and run the program. This is some explanation for learning porpouse to be read by my supervisor


# How this project has been planned

In the beginning, I used a jupyter notebook and the same couple of VCF files to develop the core of the program. Then, I run the script with different files because I expected mistakes due to the small differences between the different VCF files we generated in the lab. To solve these errors I was introducing incremental changes. Finally, when the program worked correctly, I optimized the software. For example, to select the body of the VCF files, originally I created new files with the body to load the data into a data frame. This was inefficient and very time-consuming. I improved this by selecting and importing directly the body of the files into the data frames.

Then, I adapt my code from the jupyter notebook to a script to run directly the program from the command line (these differences are explained with comments in the jupyter notebook). 

Then, immediately after I spent some time working with documentation to do not forget any important details to mention. I was adding comments while writing the code but at this point, I focused on the Readme file. 

Then, I concentrated on testing. I was testing every step while writing the code but at this point, I implement proper testing methods. This has been done in a separate folder. 



# How the env has been created

### To create the env

```
(base) monkiky@Monkikys-MacBook-Pro VCF-matcher % conda create --name vcf_matcher  python --no-default-packages

Collecting package metadata (current_repodata.json): done
Solving environment: done

## Package Plan ##

  environment location: /Users/monkiky/opt/anaconda3/envs/vcf_matcher

  added / updated specs:
    - python

The following packages will be downloaded:

    package                    |            build
    ---------------------------|-----------------
    ca-certificates-2021.7.5   |       hecd8cb5_1         113 KB
    certifi-2021.5.30          |   py39hecd8cb5_0         138 KB
    openssl-1.1.1l             |       h9ed2024_0         2.2 MB
    pip-21.2.4                 |   py37hecd8cb5_0         1.8 MB
    python-3.9.6               |       h88f2d9e_1         9.7 MB
    sqlite-3.36.0              |       hce871da_0         1.1 MB
    tzdata-2021a               |       h5d7bf9c_0         111 KB
    wheel-0.37.0               |     pyhd3eb1b0_0          32 KB
    ------------------------------------------------------------
                                           Total:        15.2 MB

The following NEW packages will be INSTALLED:

  ca-certificates    pkgs/main/osx-64::ca-certificates-2021.7.5-hecd8cb5_1
  certifi            pkgs/main/osx-64::certifi-2021.5.30-py39hecd8cb5_0
  libcxx             pkgs/main/osx-64::libcxx-10.0.0-1
  libffi             pkgs/main/osx-64::libffi-3.3-hb1e8313_2
  ncurses            pkgs/main/osx-64::ncurses-6.2-h0a44026_1
  openssl            pkgs/main/osx-64::openssl-1.1.1l-h9ed2024_0
  pip                pkgs/main/osx-64::pip-21.2.4-py37hecd8cb5_0
  python             pkgs/main/osx-64::python-3.9.6-h88f2d9e_1
  readline           pkgs/main/osx-64::readline-8.1-h9ed2024_0
  setuptools         pkgs/main/osx-64::setuptools-52.0.0-py39hecd8cb5_0
  sqlite             pkgs/main/osx-64::sqlite-3.36.0-hce871da_0
  tk                 pkgs/main/osx-64::tk-8.6.10-hb0a8c7a_0
  tzdata             pkgs/main/noarch::tzdata-2021a-h5d7bf9c_0
  wheel              pkgs/main/noarch::wheel-0.37.0-pyhd3eb1b0_0
  xz                 pkgs/main/osx-64::xz-5.2.5-h1de35cc_0
  zlib               pkgs/main/osx-64::zlib-1.2.11-h1de35cc_3


Proceed ([y]/n)? y


Downloading and Extracting Packages
openssl-1.1.1l       | 2.2 MB    | ########################################################################################## | 100% 
certifi-2021.5.30    | 138 KB    | ########################################################################################## | 100% 
tzdata-2021a         | 111 KB    | ########################################################################################## | 100% 
sqlite-3.36.0        | 1.1 MB    | ########################################################################################## | 100% 
python-3.9.6         | 9.7 MB    | ########################################################################################## | 100% 
wheel-0.37.0         | 32 KB     | ########################################################################################## | 100% 
pip-21.2.4           | 1.8 MB    | ########################################################################################## | 100% 
ca-certificates-2021 | 113 KB    | ########################################################################################## | 100% 
Preparing transaction: done
Verifying transaction: done
Executing transaction: done
#
# To activate this environment, use
#
#     $ conda activate vcf_matcher
#
# To deactivate an active environment, use
#
#     $ conda deactivate

```

#### To activate the env

```
(base) monkiky@Monkikys-MacBook-Pro VCF-matcher % conda activate vcf_matcher
(vcf_matcher) monkiky@Monkikys-MacBook-Pro VCF-matcher % pip freeze                                                   
certifi==2021.5.30
(vcf_matcher) monkiky@Monkikys-MacBook-Pro VCF-matcher % ls
Images_slides_and_stuff_to_explain_the_application	Samples
LICENSE							Test
README.md						app
Requirements.txt					test_installation.py
```

#### To install the libraries we need to run the script

```
(vcf_matcher) monkiky@Monkikys-MacBook-Pro VCF-matcher % pip install pandas
Collecting pandas
  Downloading pandas-1.3.2-cp39-cp39-macosx_10_9_x86_64.whl (11.6 MB)
     |████████████████████████████████| 11.6 MB 12.9 MB/s 
Collecting python-dateutil>=2.7.3
  Using cached python_dateutil-2.8.2-py2.py3-none-any.whl (247 kB)
Collecting pytz>=2017.3
  Using cached pytz-2021.1-py2.py3-none-any.whl (510 kB)
Collecting numpy>=1.17.3
  Downloading numpy-1.21.2-cp39-cp39-macosx_10_9_x86_64.whl (17.0 MB)
     |████████████████████████████████| 17.0 MB 14.0 MB/s 
Collecting six>=1.5
  Using cached six-1.16.0-py2.py3-none-any.whl (11 kB)
Installing collected packages: six, pytz, python-dateutil, numpy, pandas
Successfully installed numpy-1.21.2 pandas-1.3.2 python-dateutil-2.8.2 pytz-2021.1 six-1.16.0
(vcf_matcher) monkiky@Monkikys-MacBook-Pro VCF-matcher % pip install argparse
Collecting argparse
  Downloading argparse-1.4.0-py2.py3-none-any.whl (23 kB)
Installing collected packages: argparse
Successfully installed argparse-1.4.0
(vcf_matcher) monkiky@Monkikys-MacBook-Pro VCF-matcher % pip freeze          
certifi==2021.5.30
numpy==1.21.2
pandas==1.3.2
python-dateutil==2.8.2
pytz==2021.1
six==1.16.0
```

#### To see that the script work with only this libraries 
```
(vcf_matcher) monkiky@Monkikys-MacBook-Pro VCF-matcher % cd app 
(vcf_matcher) monkiky@Monkikys-MacBook-Pro app % ls
Development.ipynb	assets			run.py
(vcf_matcher) monkiky@Monkikys-MacBook-Pro app % python run.py 
File W2013397_S6.vcf contains one sample only.
File W2103016_S15.vcf contains one sample only.

 _____________________________  REPORT  ________________________________________ 

vcf 1: W2013397_S6.vcf  AND its sample name: W2013397  
vcf 2:  W2103016_S15.vcf  AND its sample name: W2103016

                                                  Homozigous: 30 
Number of positions with the same genotype: 55 
                                                  Heterozigous: 25 
                                                
                                                  
Number of positions with different genotype: 5 
                                                  


Total positions compared: 60
Percentage in common: 55/60= 0.9166666666666666

 ____________________________ END REPORT  ______________________________________
```
#### To create the  requeritments.txt with the libraries needed

```
(vcf_matcher) monkiky@Monkikys-MacBook-Pro app % ls
Development.ipynb	assets			run.py
(vcf_matcher) monkiky@Monkikys-MacBook-Pro app % cd ..
(vcf_matcher) monkiky@Monkikys-MacBook-Pro VCF-matcher % ls
Images_slides_and_stuff_to_explain_the_application	Samples
LICENSE							Test
README.md						app
Requirements.txt					test_installation.py
(vcf_matcher) monkiky@Monkikys-MacBook-Pro VCF-matcher % pip freeze > Requirements.txt
(vcf_matcher) monkiky@Monkikys-MacBook-Pro VCF-matcher % less Requirements.txt 

certifi==2021.5.30
numpy==1.21.2
pandas==1.3.2
python-dateutil==2.8.2
pytz==2021.1
six==1.16.0

```