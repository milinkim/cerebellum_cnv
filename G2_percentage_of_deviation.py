#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  7 10:07:01 2025


"""
import pandas as pd
from scipy import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.multitest import multipletests
import numpy as np
column_names = [
       'Corpus.Medullare', 'Left.Crus.I', 'Left.Crus.II', 'Left.I.III',
       'Left.IV', 'Left.IX', 'Left.V', 'Left.VI', 'Left.VIIB',
       'Left.VIIIA', 'Left.VIIIB', 'Left.X', 'Right.Crus.I',
       'Right.Crus.II', 'Right.I.III', 'Right.IX', 'Right.V',
       'Right.VI', 'Right.VIIB', 'Right.VIIIA', 'Righ.VIIIB',
       'Right.X', 'Rigt.IV', 'Vermis.IX', 'Vermis.VI', 'Vermis.VII',
       'Vermis.VIII', 'Vermis.X'] 
data=pd.read_csv('.csv')
# Calculate how many values are less than 1.96 and divide by the number of columns for each row
data['count_less_than_-1.96'] = (data[column_names] < -1.96).sum(axis=1)
data['fraction_less_than_-1.96'] = data['count_less_than_-1.96'] / len(column_names)

#filter_values = ['HC','15q11.2del', '15q11.2dup']

# Filter the DataFrame
#df = df_all[df_all['Pathogenic_CNVs'].isin(filter_values)]


# Choose specific CNV status groups to compare in the analysis
selected_cnv_statuses = ['HC', '15q11.2del', '15q11.2dup']  # Replace with the actual CNV statuses you want to compare
cnv_name='15q11.2'##########change

    
rois =  ['fraction_less_than_-1.96']  # Replace with your actual ROI column names


def anova():
    # List to collect results
    results_list = []
    
    for roi in rois:
        # Prepare data for selected Pathogenic CNV statuses
        groups_data = [data[data['Pathogenic_CNVs'] == status][roi].dropna().values for status in selected_cnv_statuses]
    
        combined_data = []
        group_labels = []
        
        for status, group in zip(selected_cnv_statuses, groups_data):
            combined_data.extend(group)
            group_labels.extend([status] * len(group))
        
        if all(len(group) > 0 for group in groups_data) and len(combined_data) == len(group_labels):
            # Perform Levene's test for homogeneity of variances
            levene_statistic, levene_p_value = stats.levene(*groups_data)
            
            try:
                # Conduct one-way ANOVA
                f_statistic, p_value = stats.f_oneway(*groups_data)
    
                # Calculate eta squared
                overall_mean = pd.Series(combined_data).mean()
                ss_total = sum((value - overall_mean)**2 for value in combined_data)
                ss_between = sum(len(group) * (group.mean() - overall_mean)**2 for group in groups_data)
                eta_squared = ss_between / ss_total
    
                # Tukey's HSD posthoc tests
                posthoc_result = pairwise_tukeyhsd(combined_data, group_labels).summary().data
                
                # Store results in a dictionary
                results_list.append({
                    'ROI': roi,
                    'F-statistic': f_statistic,
                    'ANOVA P-value': p_value,
                    'Eta Squared': eta_squared,
                    'Levene Statistic': levene_statistic,
                    'Levene P-value': levene_p_value,
                    'PostHoc': posthoc_result
                })
            except Exception as e:
                results_list.append({
                    'ROI': roi,
                    'Error': str(e),
                    'PostHoc': None
                })
        else:
            results_list.append({
                'ROI': roi,
                'Error': 'Insufficient data or data inconsistency',
                'PostHoc': None
            })
    
    # Convert results list to DataFrame
    results_df = pd.DataFrame(results_list)
    
    # Organize PostHoc results into separate columns
    posthoc_results_list = []
    for roi_result in results_df.itertuples():
        if roi_result.PostHoc:
            for row in roi_result.PostHoc:
                posthoc_results_list.append({
                    'ROI': roi_result.ROI,
                    'Group 1': row[0],
                    'Group 2': row[1],
                    'Meandiff': row[2],
                    'P-adj': row[3],
                    'Lower': row[4],
                    'Upper': row[5],
                    'Reject': row[6]
                })
    
    posthoc_results_df = pd.DataFrame(posthoc_results_list)
    
    # Save the results to CSV files
    results_df.drop(columns=['PostHoc']).to_csv(wdir+'anova_summary'+cnv_name+'.csv', index=False)
    #results_df.to_csv(wdir+'anova_summary'+cnv_name+'.csv', index=False)
    
    #posthoc_results_df.to_csv(wdir+'posthoc_summary'+cnv_name+'.csv', index=False)
    
    
    
    
    # Plotting Confidence Interval Plot for Post Hoc Tests
    for roi in posthoc_results_df['ROI'].unique():
        subset = posthoc_results_df[posthoc_results_df['ROI'] == roi]
        
        # Exclude first comparison: assuming the logic here is [Group1, Group2] vs [Group2, Group3], etc.
        subset = subset.iloc[1:]
        
        plt.figure(figsize=(10, 6))
        
        # Correctly distribute x positions
        x_positions = range(len(subset))
    
        for i, row in zip(x_positions, subset.iterrows()):
            plt.plot([i, i], [row[1]['Lower'], row[1]['Upper']], color='grey', lw=2)  # Error bar
            plt.scatter(i, row[1]['Meandiff'], color='red' if row[1]['Reject'] else 'blue', s=50)  # Mean difference
    
        plt.xticks(x_positions, labels=[f"{group1} vs {group2}" for group1, group2 in zip(subset['Group 1'], subset['Group 2'])], rotation=45, ha='right')
        plt.axhline(0, color='black', lw=2)
        plt.title(f'Post Hoc Test Results: Tukey HSD for {roi}')
        plt.xlabel('Group Comparisons')
        plt.ylabel('Mean Differences with Confidence Interval')
        plt.tight_layout()
        plt.show()
    return results_df, posthoc_results_df
        
results_df, posthoc_results_df=anova()    
####plot



group='15q11.2del' ##change
# Example of loading data into a DataFrame
# Replace with your actual data loading method
# df = pd.read_csv('your_data.csv')

# Assuming `column_names` is your list of column names needed to evaluate
column_names = [
       'Corpus.Medullare', 'Left.Crus.I', 'Left.Crus.II', 'Left.I.III',
       'Left.IV', 'Left.IX', 'Left.V', 'Left.VI', 'Left.VIIB',
       'Left.VIIIA', 'Left.VIIIB', 'Left.X', 'Right.Crus.I',
       'Right.Crus.II', 'Right.I.III', 'Right.IX', 'Right.V',
       'Right.VI', 'Right.VIIB', 'Right.VIIIA', 'Right.VIIIB',
       'Right.X', 'Rigt.IV', 'Vermis.IX', 'Vermis.VI', 'Vermis.VII',
       'Vermis.VIII', 'Vermis.X'] 
# Calculate how many values are less than 1.96 and divide by the number of columns for each row
data['count_less_than_-1.96'] = (data[column_names] < -1.96).sum(axis=1)
data['fraction_less_than_-1.96'] = data['count_less_than_-1.96'] / len(column_names)

# Segregate groups if necessary (assuming a 'group' column exists in your DataFrame)
hc_group = data[data['Pathogenic_CNVs'] == 'HC']['fraction_less_than_-1.96']
clinical_group = data[data['Pathogenic_CNVs'] == group]['fraction_less_than_-1.96']

# Perform t-test between hc and clinical groups
t_statistic, p_value = stats.ttest_ind(hc_group, clinical_group)

print(f"T-Statistic: {t_statistic}")
print(f"P-Value: {p_value}")
mean_diff_del = clinical_group.mean() - hc_group.mean()
n1, n2 = len(clinical_group), len(hc_group)
s1, s2 = clinical_group.std(ddof=1), hc_group.std(ddof=1)
pooled_std_del = np.sqrt(((n1 - 1) * s1**2 + (n2 - 1) * s2**2) / (n1 + n2 - 2))
cohens_d_del = mean_diff_del / pooled_std_del
print(f"Cohensd: {cohens_d_del}")

# Analyze the p-value to determine statistical significance
if p_value < 0.05:
    print("The difference is statistically significant.")
else:
    print("The difference is not statistically significant.")



