#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import seaborn as sns


# In[2]:


pullx_arr = np.loadtxt("pullx.xvg", dtype=float, comments=["#", "@"])

pullx_df = pd.DataFrame(pullx_arr, columns=["t", "x"]).set_index("t").apply(abs)

pullx_df


# In[3]:


numcont_arr = np.loadtxt("numcont_particle-solvent.xvg", dtype=float, comments=["#", "@"])

numcont_df = pd.DataFrame(numcont_arr, columns=["t", "numcont"]).set_index("t")

numcont_df = numcont_df.drop(numcont_df[numcont_df.numcont == 0].index)

numcont_df


# In[4]:


xnt = pullx_df.merge(numcont_df, on='t')


# In[5]:


fig = sns.jointplot(data=xnt, x="x", y="numcont", kind='hist', hue="numcont", 
                    joint_kws=dict(bins=40), marginal_kws=dict(bins=40) )
fig.savefig("numcont-x_jointplot.pdf")


# In[6]:


fig = sns.jointplot(data=xnt, x="x", y="numcont", kind='kde')
fig.savefig("numcont-x_jointplot_kde.pdf")


# In[ ]:




