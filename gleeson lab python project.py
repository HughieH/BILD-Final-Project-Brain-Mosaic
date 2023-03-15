#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 11:42:58 2023

@author: yangliu
"""
#quality control-demonstrating every sample having a high quality
#we want to focus on un-deleterious(non-pathogenetic), which are having low frequencies in GNOMAD_FREQ
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

df = pd.read_csv('./final_summary (gleeson lab data).txt', sep='\t')

#calculate the total number of alleles in a population
df["log_count"] = np.log10(df["REF_COUNT"] + df["ALT_COUNT"]) 

# Check the data type of the column
print(df['log_count'].dtype)

plt.figure(figsize=(100, 100))

# Set up the figure
fig, ax = plt.subplots()

data = df
# Create a violin plot
sns.violinplot(x=df['ID'], y=df["log_count"], data=data, hue='ID', split=False, scale='count', inner='stick', ax=ax)

# Set the x-axis label
ax.set_xlabel('samples')

# Set the y-axis label
ax.set_ylabel('log_count')

#add a horizontal line at y=3
plt.axhline(y=3, linestyle='dashed', color="blue")

plt.title('See deep length')

plt.xticks(rotation=90, fontsize=8, fontweight='bold')

# Show the plot
plt.show()