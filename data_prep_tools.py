'''A module containing tools used to prepare data for alpha diversity analyses. There are multiple tools, including those that calculate alpha diversity metrics, tools used to merge datasets and subset datasets, and tools used to plot scatterplots used for assessing independence'''

from skbio.diversity import alpha_diversity
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
from itertools import combinations

def merged_data(rarefied_OTUs, metadata):
    
    '''This function calculates two alpha diversity metrics, the Shannon Diversity Index and Observed Taxa (Sobs).
        Inputs:
        rarefied_OTUs - Rarefied OTU tables. Sample names are in first column, and OTU counts in other columns.  
        metadata - Metadata. Sample names in first column, other measurements in other columns.
        Returns:
        merged_data - The merged dataset I will use in my analyses.'''
  
    Shannon = alpha_diversity('shannon', rarefied_OTUs)
    Sobs = alpha_diversity('sobs', rarefied_OTUs)
    adiv = pd.concat([Shannon, Sobs], axis=1)
    adiv.columns = ['Shannon', 'Sobs']

    merged_data = pd.concat([adiv, metadata], axis=1)
    merged_data = merged_data.dropna(subset=["Shannon"])
   
    return merged_data
    

def scatterplots(datasets, plot_variables, labels, group):

    '''Create scatterplots assessing correlations between two groups within a particular metric and perform linear regressions. If two variables are significantly correlated, they may not be independent.
        Inputs: 
        datasets - datasets used in this analysis
        plot_variables - variables plotted
        labels - used to label datasets and plots
        group - variable that the data will be grouped by
        Returns:
        results'''

    results = []
    
    for dataset_name, dataset in datasets.items():

        for var in plot_variables:

                dataset_noNaN = dataset[[f'{var} AM', f'{var} PM']].dropna()

                x = dataset_noNaN[f'{var} AM']
                y = dataset_noNaN[f'{var} PM']

                result = stats.linregress(x, y)

                print(f'{dataset_name} {labels[var]}')
                print('p-value =', result.pvalue)

                results.append({
                    'Dataset': dataset_name,
                    'Metric': labels[var],
                    'Slope': result.slope,
                    'Intercept': result.intercept,
                    'R2': result.rvalue**2,
                    'p_value': result.pvalue,
                    'n': len(dataset_noNaN)
                })

                plt.figure()
                plt.plot(x, y, 'o')

                if result.pvalue < 0.05:
                    plt.plot(x, result.slope*x+result.intercept)

                plt.title(f"{dataset_name} {labels[var]} AM vs PM")
                plt.xlabel('AM')
                plt.ylabel('PM')
                plt.xticks(rotation=45)
                plt.show()

    return pd.DataFrame(results)
        


    
