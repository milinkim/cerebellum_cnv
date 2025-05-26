#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 17:05:50 2025

"""




#####compare bw groups for percentage of extreme 
import pandas as pd
from statsmodels.stats.proportion import proportions_ztest


# Define the extreme deviation threshold
threshold = -1.96

# Create a new column indicating extreme deviations
df['Extreme_Deviation'] = df['Total_Cerebel_Vol'] <= threshold

# Calculate counts of extreme deviations per group
deviation_counts = df.groupby('Pathogenic_CNVs')['Extreme_Deviation'].value_counts().unstack(fill_value=0)

# Prepare for pairwise comparisons
results = []

for i in range(len(deviation_counts)):
    for j in range(i + 1, len(deviation_counts)):
        group1 = deviation_counts.index[i]
        group2 = deviation_counts.index[j]
        
        # Perform the Z-test
        count = [deviation_counts.iloc[i, 1], deviation_counts.iloc[j, 1]]  # counts of extreme deviations
        nobs = [deviation_counts.iloc[i, 0] + deviation_counts.iloc[i, 1], 
                 deviation_counts.iloc[j, 0] + deviation_counts.iloc[j,1]]  # total counts
        
        z_statistic, p_value = proportions_ztest(count, nobs, alternative='two-sided')

        # Store results
        results.append((group1, group2, z_statistic, p_value))

# Output results
for result in results:
    group1, group2, z_statistic, p_value = result
    print(f'Comparison between {group1} and {group2}:')
    print(f'Z-statistic: {z_statistic:.4f}, P-value: {p_value:.4f}')
    
    # Check for significance
    if p_value < 0.05:
        print(f"Significant difference in extreme negative deviations between {group1} and {group2}.")
    else:
        print(f"No significant difference in extreme negative deviations between {group1} and {group2}.")

    print()  # New line for better readability

########################################

import pandas as pd
from scipy.stats import chi2_contingency


df=pd.read_csv('')

# Example dataframe setup
# Assuming df is your dataframe with columns: 'Pathogenic_CNVs', 'Total_Cerebellar_Vol', etc.
df = df[df['Pathogenic_CNVs'].isin(['HC', '15q11.2del', '15q11.2dup'])]
 

# Step 2: Define a function to categorize 'Total_Cerebel_Vol'
def categorize_volume(vol):
    if vol <= -1.96:
        return 'Extreme Negative'
    elif vol >= 1.96:
        return 'Extreme Positive'
    else:
        return 'Normal'

# Apply the categorization to the 'Total_Cerebel_Vol' column
df['Deviation Category'] = df['Total_Cerebel_Vol'].apply(categorize_volume)

# Step 3: Create a contingency table
# Assuming there is a 'Group' column that distinguishes between Deletion, Non-carrier, and Duplication
contingency_table = pd.crosstab(df['Pathogenic_CNVs'], df['Deviation Category'])

# Step 4: Perform the Chi-square test of independence
chi2_stat, p_value, dof, expected = chi2_contingency(contingency_table)

# Step 5: Print results
print("Contingency Table:")
print(contingency_table)
print("\nChi-square Statistic:", chi2_stat)
print("P-value:", p_value)
print("Degrees of freedom:", dof)
print("Expected frequencies:\n", expected)

# Create a summary report
summary_report = pd.DataFrame({
    'Group': contingency_table.index,
    'Extreme Negative': contingency_table['Extreme Negative'],
    'Normal': contingency_table['Normal'],
    'Extreme Positive': contingency_table.get('Extreme Positive', 0),  # Handle missing columns
    'Total': contingency_table.sum(axis=1)
})

summary_report['Percent Ext. Neg'] = (summary_report['Extreme Negative'] / summary_report['Total']) * 100
summary_report['Percent Ext. Pos'] = (summary_report['Extreme Positive'] / summary_report['Total']) * 100

print("\nSummary Report:")
print(summary_report)


####posthoc
import pandas as pd
from scipy.stats import chi2_contingency

# Assuming you already have the contingency table
contingency_table = pd.crosstab(df['Pathogenic_CNVs'], df['Deviation Category'])

# Performing the Chi-square test
chi2_stat, p_value, dof, expected = chi2_contingency(contingency_table)

# Calculate standardized residuals
residuals = (contingency_table - expected) / expected**0.5
standardized_residuals = pd.DataFrame(residuals)

# Display standardized residuals
print("Standardized Residuals:")
print(standardized_residuals)

# Determine significant differences based on standardized residuals
significant_differences = standardized_residuals[(standardized_residuals > 2) | (standardized_residuals < -2)]
print("\nSignificant Differences (standardized residuals):")
print(significant_differences.dropna(how='all'))





