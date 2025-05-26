#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 14:16:38 2025


"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np

figure_dir='cnv_figures/'

data=pd.read_csv('')
# Example DataFrame
# Specify the groups you are interested in

custom_order = ['15q11.2del', 'HC', '15q11.2dup']

# Custom colors for each group
custom_colors = {'15q11.2del': '#1f77b4', 'HC': '#ff7f0e', '15q11.2dup': '#2ca02c'}  # Assign specific colors

specific_groups = ['15q11.2del', 'HC', '15q11.2dup']

# Filter data to include only the specified groups
filtered_data = data[data['Pathogenic_CNVs'].isin(specific_groups)]

# Define lists to store the results for means and standard deviations
# Prepare lists for means and confidence intervals
# Prepare lists for means and confidence intervals
mean_values = []
confidence_intervals = []

# Z-value for 95% confidence level
z = 1.96

# Calculating means and confidence intervals for the specified groups
for pc_group in custom_order:
    group_data = data[data['Pathogenic_CNVs'] == pc_group]['Total_Cerebel_Vol']
    mean = group_data.mean()
    std = group_data.std()
    n = group_data.count()
    ci_error = z * (std / np.sqrt(n))
    mean_values.append(mean)
    confidence_intervals.append(ci_error)

# Create DataFrame for plotting
mean_ci_df = pd.DataFrame({'Pathogenic_CNVs': custom_order, 'mean': mean_values, 'ci': confidence_intervals})

# Plot using seaborn with the rocket palette
sns.set(style='white')
plt.figure(figsize=(10, 6))

# Obtain the rocket palette and swap the colors for the second and third bars
rocket_palette = sns.color_palette("rocket", len(custom_order))
custom_palette = [rocket_palette[1], rocket_palette[2], rocket_palette[0]]  # Switch second and third colors

barplot = sns.barplot(
    x='Pathogenic_CNVs', y='mean', data=mean_ci_df,
    order=custom_order, width=0.5, capsize=0.1, palette=custom_palette
)

# Define the error bar color distinct from the palette
error_bar_color = '#333333'  # Example: dark grey for error bars

# Loop through bars to add the calculated confidence intervals
for index, bar in enumerate(barplot.patches):
    bar_height = bar.get_height()
    plt.errorbar(bar.get_x() + bar.get_width() / 2, bar_height,
                 yerr=confidence_intervals[index], ecolor=error_bar_color, capsize=5)

plt.title('Mean Total Cerebel Vol with 95% Confidence Intervals by Specific Groups')
plt.xlabel('Pathogenic CNVs')
plt.ylabel('Mean Total Cerebel Vol')
y_axis_range = (-0.5, 0.1)
plt.ylim(y_axis_range)
plt.grid(False)  # Disable grid behind the bars
plt.savefig(os.path.join(figure_dir,'15q11_2_mean_diff.png'),dpi=300, bbox_inches='tight')

# Show the plot
plt.show()

