import csv
import gzip
import os
import scipy.io
import pandas as pd
import umap
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


'''
Code taken from 10x genomic: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/matrices
'''

# define MEX directory
matrix_dir = '/home/users/jvs15/MachineLearningProstate/data/GSE168668_RAW'
# read in MEX format matrix as table
mat = scipy.io.mmread(os.path.join(matrix_dir, "GSM5155457_LNCaP-RESA_matrix.mtx.gz"))
 
# list of transcript ids, e.g. 'ENSG00000243485'
features_path = os.path.join(matrix_dir, "GSM5155457_LNCaP-RESA_features.tsv.gz")
feature_ids = [row[0] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
 
# list of gene names, e.g. 'MIR1302-2HG'
gene_names = [row[1] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
 
# list of feature_types, e.g. 'Gene Expression'
feature_types = [row[2] for row in csv.reader(gzip.open(features_path, mode="rt"), delimiter="\t")]
barcodes_path = os.path.join(matrix_dir, "GSM5155457_LNCaP-RESA_barcodes.tsv.gz")
barcodes = [row[0] for row in csv.reader(gzip.open(barcodes_path, mode="rt"), delimiter="\t")]

# transform table to pandas dataframe and label rows and columns
matrix = pd.DataFrame.sparse.from_spmatrix(mat)
matrix.columns = barcodes
#matrix.insert(loc=0, column="feature_id", value=feature_ids)
matrix.insert(loc=0, column="gene", value=gene_names)
#matrix.insert(loc=0, column="feature_type", value=feature_types)

matrix = matrix.T
#Make gene names into column names and drop gene row from matrix
matrix.columns = gene_names
matrix = matrix.drop(['gene'])


# display matrix
print(matrix)
# save the table as a CSV (note the CSV will be a very large file)
#matrix.to_csv("mex_matrix.csv", index=False)

#mapper = umap.UMAP().fit()

