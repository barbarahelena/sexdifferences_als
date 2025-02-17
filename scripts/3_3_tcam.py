# TCAM pathway data
# Barbara Verhaar, barbara.verhaar@dkfz-heidelberg.de

from importlib import reload

# mprod imports
from mprod.dimensionality_reduction import TCAM
from mprod import table2tensor

# Standard libraries
import random
import datetime
import dateutil
from itertools import combinations, product, islice
from multiprocessing import Pool

# Scientific computing
import numpy as np
import pandas as pd
import re
import scipy
from sklearn.metrics import pairwise_distances, confusion_matrix, auc
from sklearn.pipeline import Pipeline
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.model_selection import StratifiedKFold
from sklearn.inspection import permutation_importance
from sklearn.preprocessing import RobustScaler, StandardScaler, MinMaxScaler, FunctionTransformer
from sklearn.ensemble import (
    AdaBoostClassifier, RandomForestClassifier, BaggingClassifier, GradientBoostingClassifier
)
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import LogisticRegressionCV, RidgeClassifier
from sklearn.metrics import plot_roc_curve

# Visualization
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# Set Seaborn style
sns.set_style('ticks')

# Permanova functions
def pairwise_melt_meta(meta, data):

    _dm = pairwise_distances(data.loc[meta['SampleID']])
    for i in range(_dm.shape[0]):
        _dm[i,i:] = None
    _dm_df = pd.DataFrame(_dm, index = meta['SampleID'].rename('sample1'), columns=meta['SampleID'].rename('sample2'))

    _dm_melt = pd.melt(_dm_df.reset_index(), id_vars=['sample1'], value_name='d').dropna()

    metacols = ['rGroup','Participant', 'SampleID']
    _dm_melt2 = meta.loc[:, metacols].merge(_dm_melt, right_on = 'sample2', left_on = 'SampleID', how = 'right').drop('SampleID', axis=1)
    _dm_melt_meta = meta.loc[:, metacols].merge(_dm_melt2, right_on = 'sample1', left_on= 'SampleID', how = 'right', suffixes = ("1","2")).drop('SampleID', axis=1)
    return _dm_melt_meta

def random_combination(iterable, r):
    "Random selection from itertools.combinations(iterable, r)"
    pool = tuple(iterable)
    n = len(pool)
    indices = sorted(random.sample(range(n), r))
    return tuple(pool[i] for i in indices)

def ss_total(data):
    return (data['d']**2).sum()

def ss_within(data):
    if 'rGroup11' in data.columns:
        # Group11 Group22 are the permuted labels of samples in Group1, Group2 resp'
        return data.query('rGroup11 == rGroup22').groupby('rGroup11')['d'].apply(lambda x: (x**2).sum()).sum()
    else:
        return data.query('rGroup1 == rGroup2').groupby('rGroup1')['d'].apply(lambda x: (x**2).sum()).sum()

def ss_between(data, ss_t):
    return ss_t - ss_within(data)

def Fstat(data, ss_t = None, N = None):
    # p = 2 as there are always 2 groups in our cases
    if N is None:
        N = len(set(data['sample1'].tolist() + data['sample2'].tolist()))
    if ss_t is None:
        ss_t = ss_total(data)
    ss_w = ss_within(data)
    ss_b = ss_t - ss_w
    return ss_b / (ss_w / (N - 2))


def gen_rperm(meta, r):
    _mice_map = meta.loc[:, ['Participant', 'rGroup']].copy()
    _mice_map.drop_duplicates(inplace=True)
    _mice_map.index = _mice_map['Participant']
    gsizes = _mice_map.groupby('rGroup').size()
    g1size, g1label, g2label = gsizes[0], gsizes.index[0], gsizes.index[-1]

    
    for combo in  random_combination(combinations(_mice_map.index, g1size), r):
        _mice_map['Group_perm'] = g2label
        _mice_map.loc[combo ,'Group_perm'] = g1label
        yield _mice_map['Group_perm'].to_dict()
        
import time 

def _permfs(args):
    __data, tv, N, permdict = args
    _data = __data.copy()
    _data['rGroup11'] = _data['Participant1'].map(permdict)
    _data['rGroup22'] = _data['Participant2'].map(permdict)
    fstat = Fstat(_data, tv, N)
    return fstat


def run_permanova(meta, data, nperms = 1000):
    _dmm = pairwise_melt_meta(meta, data)

    _dmm['rGroup11'] = _dmm['rGroup1'].copy()
    _dmm['rGroup22'] = _dmm['rGroup2'].copy()

    N = len(set(_dmm['sample1'].tolist() + _dmm['sample2'].tolist()))
    total_var = ss_total(_dmm)
    fs_obs = Fstat(_dmm, total_var, N)

    
    with Pool(processes=5) as p:
        fs_perms = p.imap(_permfs, ((_dmm, total_var, N, dd) for dd in gen_rperm(meta, nperms)))
        
        fs_perms = np.array([x for x in fs_perms])

    return ((fs_perms > fs_obs).sum() ) / (nperms ) 

# Data import
file_path = "data/pathways.csv"
data_raw = pd.read_csv(file_path, index_col=[0], sep=",")
data_raw.shape # 140 samples and 154 pathways
data = data_raw.reset_index()  # Converts MultiIndex into columns
data.rename(columns={data.columns[0]: "ID"}, inplace=True)  # Rename first index level

meta_path = "data/meta_microbiome.csv"
meta = pd.read_csv(meta_path, index_col=0, sep=",")
print(meta.head())

key_path = "data/pathwaykeys.csv"
keys = pd.read_csv(key_path, index_col=0, sep=",")
print(keys.head())

# Merge dataframes
merged_df = data.merge(meta, on="ID")  # Merge on ID, which was the first index level
print(merged_df.head())
print(merged_df.columns)
print(merged_df.index)

# filepath: /omics/groups/OE0554/internal_temp/barbara/projects/als/scripts/2_6_tcam.py
# Sort
merged_df2 = merged_df.sort_values(by=['MouseID', 'Age_ints'], ascending=[True, True])
print("Sorted DataFrame:")
print(merged_df2.head())

# Reshape
cols = ['MouseID', 'Age_ints'] + [col for col in merged_df2.iloc[:, 2:154].columns]
merged_df3 = merged_df2[cols].set_index(['MouseID', 'Age_ints'])
print("Reordered DataFrame with new index:")
print(merged_df3.head())
print(merged_df3.index)

# Normalize data with a pseudocount
pseudocount = 0.001
data_normalized = merged_df3.groupby(level='MouseID')\
                    .apply(lambda x: np.log2((x + pseudocount) / (x.query("Age_ints == 6").mean() + pseudocount)))
data_normalized = data_normalized.reset_index(level=0, drop=True)
print(data_normalized.head())
print(data_normalized.index)

# Filter timepoints 16 and smaller, excluding Age_ints 8
data_filtered = data_normalized.query("Age_ints <= 16 and Age_ints != 8")
print("Filtered data:")
print(data_filtered.head())
print(data_filtered.index)

# Check for missing data and display only columns with missing values
missing_data = data_filtered.isnull().sum()
missing_data = missing_data[missing_data > 0]
print("Columns with missing data:")
print(missing_data)

# TCAM
tensor_data, mode1_map, mode3_map = table2tensor(data_filtered, missing_flag=True)
mode1_reverse_map = {val:k for k,val in mode1_map.items()}
print(tensor_data.shape)
print(mode1_map)
print(mode1_reverse_map)

# TCAM 
tcam = TCAM(n_components=None)    # will produce minimal number of components such that
                                # total explained var > 99%

transformed_data = tcam.fit_transform(tensor_data)
df_tca = pd.DataFrame(transformed_data).rename(index = mode1_reverse_map)
rounded_expvar = np.round(100*tcam.explained_variance_ratio_, 2)
df_tca.columns = [f'F{i}:{val}%' for i,val in enumerate(rounded_expvar, start = 1)]
df_plot = meta.merge(df_tca, left_on = 'MouseID', right_index=True)
f1 = df_tca.columns[0]
f2 = df_tca.columns[1]

fig,ax  = plt.subplots(1,1, figsize=[4,4], dpi = 500)
sns.scatterplot(data = df_plot, x = f1, y = f2 , hue='GenotypePerSex', ax = ax, edgecolor = 'k', linewidths=.1)
ax.legend(loc='center left', fontsize = 6, fancybox = False, framealpha = 1, edgecolor = 'k', bbox_to_anchor=(1, 0.5))
plt.savefig('results/pathways/tcam/scatterplot_paths.png', format='png', dpi=300, bbox_inches='tight')

df_plot.to_csv('results/pathways/tcam/df_plot_paths.csv', index=False)

# Permanova
n_factors = (np.cumsum(tcam.explained_variance_ratio_) < .2).sum() + 1  # PC for 20% variance
perm_res = pd.DataFrame(columns=['g1', 'g2', 'p'])
# Ensure meta_perm is correctly defined
meta_perm = meta
meta_perm['Participant'] = meta_perm['MouseID']
meta_perm['SampleID'] = meta_perm['MouseID']
print("meta_perm head:")
print(meta_perm.head())

# Check unique values in 'GenotypePerSex'
unique_values = meta_perm['GenotypePerSex'].unique()
print("Unique values in 'GenotypePerSex':")
print(unique_values)

# Check if there are enough unique values to generate combinations
for g1, g2 in combinations(unique_values, 2):
    print(g1)
    print(g2)
    meta_comp = meta_perm.loc[meta_perm['GenotypePerSex'].isin([g1, g2]) & (meta_perm['Age_ints'] == 6)].copy()
    meta_comp.rename(columns={'GenotypePerSex': 'rGroup'}, inplace=True)
    table_comp = df_tca.loc[meta_comp['MouseID']].copy()
    
    n_all = meta_comp.loc[meta_comp['rGroup'].isin([g1, g2])].shape[0]
    n_1 = meta_comp.loc[meta_comp['rGroup'].isin([g1])].shape[0]
    n_2 = meta_comp.loc[meta_comp['rGroup'].isin([g2])].shape[0]
    n_perm = min(1000, n_all)  # Ensure n_perm is not larger than the number of available samples
    
    # Ensure run_permanova is correctly defined and used
    try:
        p = run_permanova(meta_comp, table_comp.iloc[:, :n_factors].copy(), nperms=n_perm)
        print("PERMANOVA p-value for {} vs {}: {}".format(g1, g2, p))
        perm_res.loc[perm_res.shape[0], :] = [g1, g2, p]
    except Exception as e:
        print("Error running PERMANOVA for {} vs {}: {}".format(g1, g2, e))

print(perm_res)
perm_res.to_csv("results/pathways/tcam/perm_res_paths.csv")

# Loadings
loadings = tcam.mode2_loadings
df_loadings = pd.DataFrame(loadings, index =  data_normalized.iloc[:,0:].columns)
df_loadings.columns = [f'F{i}:{val}%' for i,val in enumerate(rounded_expvar, start = 1)]
# df_loadings = df_loadings.merge(key_path, on="path")
df_loadings.to_csv('results/pathways/tcam/df_loadings_paths.csv', index=True)

# PC
loadings_f1 = df_loadings.iloc[:,0].sort_values().copy()
top_bottom = [ loadings_f1.nlargest(10), loadings_f1.nsmallest(10) ]
loadings_f1 = pd.concat(top_bottom)
loadings_f1 = loadings_f1.sort_values()
loadings_f1 = loadings_f1.rename('F1')
loadings_f1 = loadings_f1.reset_index()

fig,ax  = plt.subplots(1,1, figsize=[1.2,2], dpi = 500)
ax.barh(y = loadings_f1.query('F1 < 0').index, width = loadings_f1.query('F1 < 0')['F1'] 
       , linewidth = .25, edgecolor = 'k')
ax.barh(y = loadings_f1.query('F1 > 0').index.astype(int), width = loadings_f1.query('F1 > 0')['F1']
        , linewidth = .25, edgecolor = 'k')
ax.set_yticks([-.7,loadings_f1.index.max() + 3])
ax.set_yticklabels([])

loadings_f1.index.names
for i, row in loadings_f1.iterrows():
    sname = str(row['index'])
    sname = sname.replace(";__", "")
    sname = sname.replace(";s__", "_")
    sname = re.sub(";s__$", "", sname)
    sname = re.sub('\w__[;]*','', sname)
    sname = re.sub("[;_]*$","",  sname)
    sname = sname.split(";")[-1]
    ax.text(x=0.7, y=i, s = sname, fontdict={'size':3})
    
ax.yaxis.set_tick_params(size = 3, length = 400,direction = 'inout')
sns.despine( top=True,bottom=True,right=True,trim=True)
ax.spines['left'].set_position('zero')
plt.savefig('results/pathways/tcam/barplot_pc1_paths.png', format='png', dpi=300, bbox_inches='tight')

# PC 2
loadings_f2 = df_loadings.iloc[:,1].sort_values().copy()
top_bottom = [ loadings_f2.nlargest(10), loadings_f2.nsmallest(10) ]
loadings_f2 = pd.concat(top_bottom)
loadings_f2 = loadings_f2.sort_values()
loadings_f2 = loadings_f2.rename('F2')
loadings_f2 = loadings_f2.reset_index()

fig,ax  = plt.subplots(1,1, figsize=[1.2,2], dpi = 500)
ax.barh(y = loadings_f2.query('F2 < 0').index, width = loadings_f2.query('F2 < 0')['F2'] 
       , linewidth = .25, edgecolor = 'k')
ax.barh(y = loadings_f2.query('F2 > 0').index.astype(int), width = loadings_f2.query('F2 > 0')['F2']
        , linewidth = .25, edgecolor = 'k')
ax.set_yticks([-.7,loadings_f2.index.max() + 3])
ax.set_yticklabels([])

loadings_f2.index.names
for i, row in loadings_f2.iterrows():
    sname = str(row['index'])
    sname = sname.replace(";__", "")
    sname = sname.replace(";s__", "_")
    sname = re.sub(";s__$", "", sname)
    sname = re.sub('\w__[;]*','', sname)
    sname = re.sub("[;_]*$","",  sname)
    sname = sname.split(";")[-1]
    ax.text(x=0.7, y=i, s = sname, fontdict={'size':3})
    
ax.yaxis.set_tick_params(size = 3, length = 400,direction = 'inout')
sns.despine( top=True,bottom=True,right=True,trim=True)
ax.spines['left'].set_position('zero')
plt.savefig('results/pathways/tcam/barplot_pc2_paths.png', format='png', dpi=300, bbox_inches='tight')