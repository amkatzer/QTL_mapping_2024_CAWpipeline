# QTL_mapping_2024_CAWpipeline
QTL mapping of Penstemon barbatus x Penstemon neomexicanus F2 population using Carrie Wessinger's protocol starting after the concatenation steps

## Step 1: concatenate files together
Concatenate any sequencing files together from multiple plates. Each sample should have the format "NB#" without a zero before the number (e.g. NB001 is incorrect, NB1 is correct).
Script to run for F2s: cat.fastq.sh
Script to run for parents: 

## Step 2: remove any samples from concatenation script output where there is no sequencing data. 
Delete anything with no data in the file to clean up the folder.

## Step 3: make a sample list from the remaining outputs from the concatenation script
