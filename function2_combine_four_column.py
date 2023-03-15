#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 10:26:16 2023

@author: yangliu
"""

import pandas as pd
import numpy as np

# Load the Excel file
df = pd.read_csv('BILD-Final-Project-Brain-Mosaic/final_summary (gleeson lab data).txt', sep='\t')

# print(df.columns)
# df.describe()
# Subset the rows based on a condition: dimerge the graph for each samples and if each sample has an average less than three, we will need to filter.
# & exclude them

df = df[~(df['ID'] == '7669_L_sec1_Put_MSN')]

# don't think we need this? @robert
# np.linspace(0,147,10)

# Print the subset of rows
#print(exclude_df)

# Add new columns for Chromosome, Position, REF, ALT
# df["CHROM_POS_REF_ALT"] = None #overlap, several variants in same chromosome or position so that u need to find out the combination of them 

df.insert(5, 'CHROM_POS_REF_ALT', None)

# mosaic is if the variant is present in some tissues but not in others, it indicates the variant is mosaic.
# df["Mosaic"] = None #distinguish whether if it is mosaic

df.insert(6, 'Mosaic', None)

#pastes together the values of four columns and assign them to one column
df["CHROM_POS_REF_ALT"] = df['CHROM'].astype(str) + '-' + df['POS'].astype(str) + '-' + df['REF'].astype(str) + '-' + df['ALT'].astype(str)

# Check if any of the CHROM_POS_REF_ALT are identical, if so remove them from the dataframe
duplicates = df['CHROM_POS_REF_ALT'].duplicated()
print(duplicates.head) # check which rows are duplicates

df = df[~duplicates]

# Relocate the new column after the 'ALT' column
df.to_csv('final_summary_anno3_hugh.txt', sep='\t', index=False)
df.to_csv('final_summary_anno3_hugh.csv', index=False)

