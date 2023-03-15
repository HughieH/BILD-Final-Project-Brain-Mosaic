#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 19:58:19 2023

@author: yangliu
"""
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats


# Read the TSV file into a data frame
df = pd.read_csv('./final_summary_anno3.txt', sep='\t', index_col=0)

print(df.head)

df.shape
# function define

variants = df.CHR_POS_REF_ALT.unique()

#summarize maf
def condense_maf_data(dbsm_data):
    final_data = dbsm_data.loc[df.CHR_POS_REF_ALT == variants[0], ["ID", "MAF"]]
    final_data.columns = ["ID", variants[0]]
    for variant in variants[1:len(variants)]:
        data = dbsm_data.loc[df.CHR_POS_REF_ALT == variant, ["ID", "MAF"]]
        data.columns = ["ID", variant]
        final_data = pd.merge(final_data, data, on="ID", how="outer")
    final_data = final_data.set_index(final_data.ID)
    final_data = final_data[final_data.columns[1:]]
    return final_data 

#compute correlations among all variant pairs
def compute_pairwise_corr(maf_df):
    coef_matrix = np.zeros([len(variants), len(variants)])
    pval_matrix = np.zeros([len(variants), len(variants)])
    for i in range(len(variants)):
        var_i = variants[i]
        mafs_i = maf_df[var_i].values
        for j in range(len(variants)):
            var_j = variants[j]
            mafs_j = maf_df[var_j].values
            data = pd.DataFrame(np.vstack([mafs_i, mafs_j]).T).dropna()
            coef, pval = scipy.stats.pearsonr(data[0], data[1])
            coef_matrix[i,j] = coef
            pval_matrix[i,j] = pval
    coef_df = pd.DataFrame(coef_matrix, columns = variants)
    coef_df = coef_df.set_index(variants)
    pval_df = pd.DataFrame(pval_matrix, columns = variants)
    pval_df = pval_df.set_index(variants)
    return coef_df, pval_df



#summarize mafs for sorted population data and bulk data

maf_df = condense_maf_data(df)

coef_df, pval_df = compute_pairwise_corr(maf_df)

#remove NA rows and columns
coef_df = coef_df.dropna()
variants = coef_df.index
coef_df = coef_df[variants]

#define the color coding 

pf_color_dict = {True: "#5C3763", False: "#FFFFFF"}
f_color_dict = {True: "#88527B", False: "#FFFFFF"}
p_color_dict = {True: "#AF738D", False: "#FFFFFF"}
o_color_dict = {True: "#CF9DA4", False: "#FFFFFF"}
t_color_dict = {True: "#E8CCC7", False: "#FFFFFF"}

cerebellum_color_dict = {True: "#173E20", False: "#FFFFFF"}
heart_color_dict = {True: "#3D7247", False: "#FFFFFF"}
liver_color_dict = {True: "#679D72", False: "#FFFFFF"}
kidney_color_dict = {True: "#A0C7A7", False: "#FFFFFF"}

#load data

#value mapping between variants and some columns
pf_dict = dict(zip(df.CHR_POS_REF_ALT, df.IN_PF))
f_dict = dict(zip(df.CHR_POS_REF_ALT, df.IN_F))
p_dict = dict(zip(df.CHR_POS_REF_ALT, df.IN_P))
o_dict = dict(zip(df.CHR_POS_REF_ALT, df.IN_O))
t_dict = dict(zip(df.CHR_POS_REF_ALT, df.IN_T))
cerebellum_dict = dict(zip(df.CHR_POS_REF_ALT, df.SET_IN_CEREBELLUM))
heart_dict = dict(zip(df.CHR_POS_REF_ALT, df.SET_IN_HEART))
liver_dict = dict(zip(df.CHR_POS_REF_ALT, df.SET_IN_LIVER))
kidney_dict = dict(zip(df.CHR_POS_REF_ALT, df.SET_IN_KIDNEY))



'''map creating'''


#cluster map-sorted data

variants = coef_df.columns
pf_colors = variants.map(pf_dict).map(pf_color_dict)
f_colors = variants.map(f_dict).map(f_color_dict)
p_colors = variants.map(p_dict).map(p_color_dict)
o_colors = variants.map(o_dict).map(o_color_dict)
t_colors = variants.map(t_dict).map(t_color_dict)
cerebellum_colors = variants.map(cerebellum_dict).map(cerebellum_color_dict)
heart_colors = variants.map(heart_dict).map(heart_color_dict)
liver_colors = variants.map(liver_dict).map(liver_color_dict)
kidney_colors = variants.map(kidney_dict).map(kidney_color_dict)

clustermap_sorted = sns.clustermap(coef_df, figsize=(50,50), cmap="RdBu_r",
                row_colors=[cerebellum_colors, heart_colors, liver_colors, kidney_colors,
                            pf_colors, f_colors, p_colors, o_colors, t_colors],
              colors_ratio=0.005, dendrogram_ratio=0.1)
plt.savefig('population_clustermap.svg')     

