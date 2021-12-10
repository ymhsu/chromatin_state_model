## Environment and software installation
All of the analyses were performed under Linux (Ubuntu 18.04.6 LTS).

The shell scripts use bash.

To perform the analyses, one needs to have a number of utilities installed:  
"R" (we use R Studio version 1.1.463 with R version 3.6.3.)  
The necessary libraries are included therein  
"bedtools" (we use version of bedtools v2.26.0, free from https://github.com/arq5x/bedtools2/releases)   
"bigWigToBedGraph" (free from http://hgdownload.cse.ucsc.edu/admin/exe/)  

## Repository organization
The files are organized in 3 directories as follows:  


### data directory: 
contains names of datasets to be downloaded for each figure as well as the CO datasets of the publications Rowan et al. 2020 and Blackwell et al 2020.  


### script directory:
Contains 5 shell scripts used to download raw data and produce intermediate files needed by R scripts.  
Contains 6 R codes used for performing the analyses published in the paper (each code can be operated independently, just make sure the shell scripts are executable).  
All these scripts contain documentation.   
Note that the R codes execute internally all the required shell scripts.   
The R scripts also install required dependencies. 


### analysis_output:
The results are written here when running the scripts


## Cloning the repository
The simplest way to reproduce the analyses is to use git to clone the directories via:
```bash
git clone https://github.com/ymhsu/chromatin_state_model.git
``` 
This will produce a directory "chromatin_state_model" containing this Readme file and the 3 subdirectories described above.


## Illustrative examples
### The interactive mode
All of R codes are initially written for the R-studio environment.
After cloning the repository by "git clone", one can create a new project of R-studio from the directory "chromatin_state_model".
This will automatically set the directory "chromatin_state_model" as the working directory for running all R codes.


### The noninteractive mode
#### Reproducing Figures 1, Figure S1 and Table S2 of the paper
From the directory "chromatin_state_model/", execute:
```bash
R CMD BATCH ./script/Fig1_S1_final.R
``` 

Note: the shell script "Fig1_9features_intersect_14bins_f.sh", used in this R code, produces the required data files for 14 different bin sizes ranging from 0.5 kb to 1000 kb. You can change that both in the shell script and this R code to your favorite bin size(s).

#### Reproducing Figure 2 and S2-S6 of the paper
From the directory "chromatin_state_model/", execute: 
```bash
R CMD BATCH ./script/Fig2_new.R
```

#### Reproducing Figure 3 of the paper
From the directory "chromatin_state_model/", execute: 
```bash
R CMD BATCH ./script/Fig3_The_effect_IR_sizes.R
```

#### Reproducing Figure 4 and S7 of the paper
From the directory "chromatin_state_model/", execute: 
```bash
R CMD BATCH ./script/Fig4_new.R
```

#### Reproducing Figure 5, Figure S8, Figure S9, Table S4 and Table S5 of the paper
From the directory "chromatin_state_model/", execute: 
```bash
R CMD BATCH ./script/Fig5_new_final.R
```

#### Reproducing Table S3 of the paper
From the directory "chromatin_state_model/", execute: 
```bash
R CMD BATCH ./script/the_comparison_AIC_BIC_of_models.R
```
