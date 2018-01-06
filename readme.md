# Change Point Detection Toolbox

Change point detection algorithms for segmenting time series data.

## Overview

Segmentation performance of eight different algorithms is tested in this repo. 

## Dependencies

HRVToolbox for real data analysis part

## Configuration

+ Install R Studio to implement BCP method. Add 'bcp' package which is located here:
https://cran.r-project.org/web/packages/bcp/index.html

Add also stingr package:
install_package("stringr")

+ Install NAG MATLAB toolbox from following address:
https://www.nag.com/content/downloads-nag-toolbox-matlab-versions

+ The executable "rrgenV3" is compiled for macOS. For other operating systems, follow these steps:

1. cd into rrGenerator folder

2. Delete `rrgenV3` since this executable is for mac

3. Compile `rrgenV3.c` into executable:

    `gcc -Wall rrgenV3.c -lm -o rrgenV3`

4. Test the executable by calling `rrgenV3` from the terminal, using the number `1` as a seed, and setting max elements to `200`:  
    `./rrgenV3 1 200` % For mac computer  <br />
    `:rrgenV3 1 200` % For windows computer <br />

You should see output printed to the terminal. <br />

## Algorithms

+ **Modified Bayesian Online Changepoint Detection (Modified BOCD):** [`bocpd_revised.m`]() 
+ **Original Bayesian Online Changepoint Detection (Original BOCD):** [`bocpd_original.m`]() 
+ **Bayesian Blocks (BBLOCKS):** [`bblocks.m`]()
+ **Pruned Exact Linear Time (PELT1)** from NAG MATLAB toolbox 
+ **Matlab function (PELT2) 'findchangepts'**
+ **Binary Segmentation (BiS)** from NAG MATLAB toolbox 
+ **Recursive Mean Difference Maximization (RMDM)** : [`rmdm.m`]()
Includes 2 versions, [`rmdm.m_v1`]() is written by Maxim Osipov and [`rmdm.m`]() is the 
modified version.
+ **Bayesian Analysis of Changepoint Problems (BCP)** : [`bcp_testRRGEN.R`]()

## Real data analysis

Physionet capslpdb database is located here:
https://www.physionet.org/pn6/capslpdb/
+ Download hypnograms from above address. 
+ Download EKG files in csv format.  <br />
ex: rdsamp -r capslpdb/n5.edf -H -f 0 -ps -c >n5.csv  
 
+ Arrange files in same format as CPD/Toolbox/testData/physionet_capslpdb. EKG file and hypnogram
of same subject should be named the same. ex.: n1.txt 

## Possible errors & how to troubleshoot them

+ /bin/bash: Rscript: command not found <br />
Add where R is located to your MATLAB environment <br />
setenv('PATH', [getenv('PATH'),':','/usr/local/bin']); 

+ In artificial data analysis, no data is generated
Make sure your rrgen_2003_v2 executable is suitable for your computer type. 
If not, delete this and re compile as decribed above.

## Useful Links

Related to setting up parameters of NAG functions: <br />
PELT: http://www.nag.com/numeric/mb/nagdoc_mb/manual_25_1/html/g13/g13naf.html <br />
BiS: http://www.nag.com/numeric/mb/nagdoc_mb/manual_25_1/html/g13/g13ndf.html <br />
