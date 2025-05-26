#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 10:57:20 2025

"""

import pandas as pd
import numpy as np
import statsmodels.api as sm



#plot deviation of voxels for all groups

import random 
import glob
import pandas
import os
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import joypy
from sklearn.model_selection import train_test_split
from pcntoolkit.normative import estimate, evaluate
from pcntoolkit.util.utils import create_bspline_basis, compute_MSLL
import pcntoolkit as pcn
import pickle

import SUITPy.flatmap as flatmap
import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np

from math import sqrt
import pingouin as pg
from scipy import stats
from scipy.stats import ttest_ind
from colorspacious import cspace_converter

wdir="/"
#load mask and example image and test
mask_image='t1_bin.nii' #mask thresholded at 0.2
example_image='t1_n4_mni_seg_post_inverse.nii'
cov_dir=''
figure_dir=''
te_cnv=pd.read_csv(')#select right df_te
te=pd.read_csv(')
total=pd.read_csv(/Z_estimate.txt', header=None)

te=pd.merge(total, te, left_index=True, right_index=True)

#qc
qc=pd.read_csv('/')
cnv_names=qc['Pathogenic_CNVs'].value_counts()

dose_15q112=['15q11.2del','15q11.2dup']


te_cnv=pd.read_csv(')
te=pd.read_csv(') 

df_te=pd.read_csv(cov_dir+'cnv_aca_cov_te_clinical.csv')
df_te=df_te[['ID']]
#merge_tr=pd.concat([bl_cov_tr,merge_df_te])

total.columns = ['Total_Cerebel_Vol']

df_te=pd.merge(total, df_te, left_index=True, right_index=True)
te_cnv['ID']=te_cnv['ID'].astype(str)
te_cnv=te_cnv.merge(df_te, how= "inner", on =["ID"])
te=te.merge(df_te, how= "inner", on =["ID"])

columns_to_plot= [
       'Corpus.Medullare', 'Left.Crus.I', 'Left.Crus.II', 'Left.I.III',
       'Left.IV', 'Left.IX', 'Left.V', 'Left.VI', 'Left.VIIB',
       'Left.VIIIA', 'Left.VIIIB', 'Left.X', 'Right.Crus.I',
       'Right.Crus.II', 'Right.I.III', 'Right.IX', 'Right.V',
       'Right.VI', 'Right.VIIB', 'Right.VIIIA', 'Right.VIIIB',
       'Right.X', 'Rigt.IV', 'Vermis.IX', 'Vermis.VI', 'Vermis.VII',
       'Vermis.VIII', 'Vermis.X', 'Total_Cerebel_Vol'] 



#cnv_names=['15q11.2del','15q11.2dup']


def clinical_cohorts( te_cnv, c_name):
    cnv_df=te_cnv[te_cnv['Pathogenic_CNVs'].isin(c_name)]
  
    return cnv_qc


def hc_cohort( te):
  
    te = te[te['site'].str.startswith('ukb')] #only get ukb-so noncarriers
    return te


def determine_copy_number_status(value):
    if 'del' in value:
        return 1
    elif 'dup' in value:
        return 3
    else:
        return None  # or some other default, like 0 or -1

def run_dose(df,cnv_dose_name):
    results = []
    for column in columns_to_plot: #
        if column in columns_to_plot:
            alpha = 0.05
            num_rois = (19)
            num_tests_per_roi = 2  # deletion vs non-carrier and duplication vs non-carrier
            total_tests = num_rois * num_tests_per_roi
            bonferroni_alpha = alpha / total_tests

            X = df['Copy_Number']
            y = df[column]
            
            # Add a constant term for the intercept
            X = sm.add_constant(X)
            
            # Fit a linear regression model
            #model = sm.OLS(y, X).fit()
            model = sm.OLS(y, X).fit(cov_type='HC3')  # HC3 for heteroskedasticity-consistent SEs
            #print("Linear Regression Results:")
            #print(model.summary())
            
            #pairwise T-Tests with Welch's correction
            deletion = df[df['Copy_Number'] == 1][column]
            non_carrier = df[df['Copy_Number'] == 2][column]
            t_stat_del, p_val_del = stats.ttest_ind(deletion, non_carrier, equal_var=False)

            # Duplication vs Non-carrier
            duplication = df[df['Copy_Number'] == 3][column]
            t_stat_dup, p_val_dup = stats.ttest_ind(duplication, non_carrier, equal_var=False)
            
            t_stat_dup_del, p_val_dup_del = stats.ttest_ind(duplication, deletion, equal_var=False)

            sig_del_nc = p_val_del < bonferroni_alpha
            sig_dup_nc = p_val_dup < bonferroni_alpha
            sig_del_dup = p_val_dup_del < bonferroni_alpha


            #plot
            plt.gcf().set_dpi(300)
            # Boxplot with jittered datapoints
            plt.figure(figsize=(10, 6))
    
            sns.boxplot(x='Copy_Number', y=column, data=df, palette="rocket", linewidth=2.5, fliersize=3)
            sns.stripplot(x='Copy_Number', y=column, data=df, color='black', jitter=True, alpha=0.1, size=1)
            
            #asterick
            # Asterisks for significance
            # Asterisks and brackets for significance
            x_positions = [0, 1, 2]  # ['deletion', 'non-carrier', 'duplication']
            y_max = df[column].max() + 0.5
            y_step = 0.2  # height increment for significance lines
            
            # Plot significance asterisks and brackets
            def plot_significance_brackets(x1, x2, y, significance, ax):
                ax.plot([x1, x1, x2, x2], [y, y+y_step, y+y_step, y], lw=1.5, color='black')
                if significance:
                    ax.text((x1+x2)*.5, y+y_step, "*", ha='center', va='bottom', color='black', fontsize=15)
                return y + y_step + 0.1
            
            ax = plt.gca()
            if sig_del_nc:
                y_max = plot_significance_brackets(x_positions[0], x_positions[1], y_max, sig_del_nc, ax)
            
            if sig_dup_nc:
                y_max = plot_significance_brackets(x_positions[1], x_positions[2], y_max, sig_dup_nc, ax)
            
            if sig_del_dup:
                y_max = plot_significance_brackets(x_positions[0], x_positions[2], y_max, sig_del_dup, ax)
            
                        #sns.boxplot(x='copy_number_status', y=column, data=df, showfliers=False)
            #sns.stripplot(x='copy_number_status', y=column, data=df, color='black', jitter=True, alpha=0.5)
            plt.grid(visible=True, color='lightgray', linestyle='--', linewidth=0.5)
           
            plt.title(f'{cnv_dose_name} distribution of {column}')
            plt.xlabel('CNV')
            plt.ylabel('Z-score')
            plt.ylim(df[column].min() - 0.5, y_max + 1)  # Extend y_max to ensure the asterisks aren't clipped
            plt.xticks(ticks=[0, 1, 2], labels=['deletion', 'non-carrier', 'duplication'])

            plt.tight_layout()

            #plt.show()
            plt.savefig(os.path.join(figure_dir,cnv_dose_name+'_'+column +'_'+'_boxplot.png'))

            # Print the summary of the regression model
            #print(model.summary())
   
            params = model.params[0]  # Coefficients
            pvalues = model.pvalues[1]  # p-values
            rsquared = model.rsquared 
            results.append({
                'Column': column,
                'param': params,
                'P-value': pvalues,
                'rsquared': rsquared,
                'T-stat Del vs Non-carrier': t_stat_del,
                'P-value Del vs Non-carrier': p_val_del,
                'Bonferroni Significant Del': p_val_del < bonferroni_alpha,
                'T-stat Dup vs Non-carrier': t_stat_dup,
                'P-value Dup vs Non-carrier': p_val_dup,
                'Bonferroni Significant Dup': p_val_dup < bonferroni_alpha,
                'T-stat Dup vs Del': t_stat_dup_del,
                'P-value Dup vs Del': p_val_dup_del,
                'Bonferroni Significant Dup vs Del': p_val_dup_del < bonferroni_alpha
            })
    # Convert the list of results into a DataFrame
    results_df = pd.DataFrame(results)
    return results_df

df1 = clinical_cohorts(te_cnv, dose_15q112)##change#########
df2 = hc_cohort(te)
df2['Copy_Number']=2

df1['Copy_Number'] = df1['Pathogenic_CNVs'].apply(determine_copy_number_status)

df=pd.concat([df1, df2])

results_df=run_dose(df,cnv_dose_name)
#results_df.to_csv(wdir+'/dose_effect_'+cnv_dose_name+'.csv')
    

