import sys
import os
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
import pandas as pd
import seaborn as sns
import numpy as np
import scipy.stats as st
from sklearn.cluster import AgglomerativeClustering

def findCutoffs(df):
    d_points=df['TMB-Total-Variants'].values.reshape([-1,1])
    t_points=np.reshape(range(0,350),(-1,1))

    model = AgglomerativeClustering(n_clusters=3).fit(d_points)
    df['hclust'] = model.fit_predict(d_points)

    centers={}
    centers[0]=df[df['hclust']==0]['TMB-Total-Variants'].mean()
    centers[1]=df[df['hclust']==1]['TMB-Total-Variants'].mean()
    centers[2]=df[df['hclust']==2]['TMB-Total-Variants'].mean()

    predict =[]
    for point in t_points:
        predict.append(np.argmin([np.abs(point-centers[0])[0],np.abs(point-centers[1])[0],np.abs(point-centers[2])[0]]))

    cutoffs =[]
    start = predict[0]
    for i in range(len(predict )):
        if predict[i]!=start:
            cutoffs.append(i+1)
            start = predict[i]

    return cutoffs
df = pd.read_csv("tmb_validation_samples_unique.csv")
cutoffs = []
for i in range(1500):
    temp = findCutoffs(df.sample(n = df.shape[0], replace = True))
    if len(temp)>1:cutoffs.append(temp)

df_cut = pd.DataFrame(cutoffs)
df_cut2 = df_cut.sample(n=1000,replace = False)
df_cut2.to_csv("temp_cutoff.csv")
