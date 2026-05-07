This is Lindsay Yue’s final project for MS 263 Spring 2026.

Repository Layout:
This repository contains an empty folder called Figures, where .png files of all figures produced will go, the initial datasets used in .csv or .tsv format, the main project notebook in an .ipynb file, and .py files with functions used in the project. 

Data information:
All data is from experiments performed by the Jue lab at CSUMB in 2019 and 2024. Raw sequence data  was  processed using QIIME2 bioinformatic pipelines, and one output was the OTUs present in each sample. In R, a phyloseq object was created using QIIME2 outputs, and OTU tables were created. More about this data can be learned by contacting Lindsay Yue (lindsay.yue@sjsu.edu) or Dr. Nathaniel Jue (njue@csumb.edu). For this project, OTU tables and sampling experiment metadata were used, and can be found in this repository in the Data folder. 

Import the data as a pandas dataframe, and then proceed with the code in Analysis.ipynb.
