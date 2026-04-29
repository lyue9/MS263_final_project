'''A module for calculating the two alpha diversity metrics, Shannon diversity and Sobs, and connecting that dataset with the metadata and rarefied OTU table'''

from skbio.diversity import alpha_diversity
import pandas as pd

def final_dataset(rarefied_OTUs, metadata):
    Shannon = alpha_diversity('shannon', rarefied_OTUs)
    Sobs = alpha_diversity('sobs', rarefied_OTUs)
    adiv = pd.concat([Shannon, Sobs], axis=1)
    adiv.columns = ['Shannon', 'Sobs']

    final_data = pd.concat([adiv, metadata], axis=1)
    final_data.dropna(subset=["Shannon"])
   
    return final_data

    '''Input:
       rarefied_OTUs - Rarefied OTU tables. Sample names are in first column, and OTU counts in other columns.  
       metadata - Metadata. Sample names in first column, other measurements in other columns.
       Returns:
       final_data - The final dataset I will use in my analyses.'''
  
