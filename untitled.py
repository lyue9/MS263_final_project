'''A module used making violin plots showing the Shannon diversity and observed taxa from my water and sediment samples grouped various ways, grouped by year and grouped by port and dosing'''

import matplotlib.pyplot as plt
import seaborn as sns

def violin_plots_year(datasets):
    Shannon = alpha_diversity('shannon', rarefied_OTUs)
    Sobs = alpha_diversity('sobs', rarefied_OTUs)
    adiv = pd.concat([Shannon, Sobs], axis=1)
    adiv.columns = ['Shannon', 'Sobs']

    merged_data = pd.concat([adiv, metadata], axis=1)
    merged_data = merged_data.dropna(subset=["Shannon"])
   
    return merged_data

    '''Input:
       rarefied_OTUs - Rarefied OTU tables. Sample names are in first column, and OTU counts in other columns.  
       metadata - Metadata. Sample names in first column, other measurements in other columns.
       Returns:
       merged_data - The final dataset I will use in my analyses.'''
  