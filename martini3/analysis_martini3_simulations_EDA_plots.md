---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.5.2
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

```python
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#from cycler import cycler
#from sklearn.decomposition import PCA
from sklearn import preprocessing

#from IPython.display import display
#import MDAnalysis as mda
%pylab inline
```

# Read-in the dataset from previous step

Here basic load-in and table head

```python
# read in the csv file
data = pd.read_csv("dataFrame_all_sims.csv")
data.head()
```

```python jupyter={"outputs_hidden": false} run_control={"marked": true}
# Number of rows
print("Data have {} records".format(data.shape[0]))

# Column names
print("\nColumn names are\n{}\n".format(data.columns.tolist()))

# Data types
print("Data types of the columns are:\n {}".format(data.dtypes))
```

More detailed information about the dataset:

```python
data.info()
```

There seem to be some missing values in Compressibility column,
and there is an extra column "Unnam:0", 
which is a duplicate index. 

Let's get rid of it -- turning it into "ID"

```python
data.rename(columns={"Unnamed: 0" : "ID"}, inplace=True)
data.stack()
data.info()
```

And Now I turn "ID" column into a (simple!) Multiindex:

```python
# create multiindex from columns 1-8
midata = data.set_index(list(data.columns[:1]))
midata.head()
```

```python
# confirming that the multiindex is formed as intended:
print(midata.index.names)
```

## Missing values

Now, 
there are several missing values in Compressibility column.

Track them down and decide on what to do about them.
(maybe nothing?)

```python
# Locate the indices of the records with missing values
indices_of_missing_compr_data = midata.isna().query('compress == True').index
# printout the slice of the midata df with the missing data
midata.loc[indices_of_missing_compr_data]
```

Ha, the data arise from a separate "confirming" simulation. 
There already is an equivalent set of these simulations including such data:

```python
# using the slicing based on the above data set to
# provide a list of simulations with only particle "S"
# and with saturation index = 0.0
midata.loc[midata.particle == "S"].query("`satur index` == 0.0")
```

So, 
the simplest thing to do with such a data set is to skip those lines
as the values quite agree together. 

However, it would be even better to make an average of the duplicate records...
(let's keep that for future)

For now, I will just rename the weirdly call values (appended with "(1)")
with their "normal" equivalents 
and will use tham alongside the other values. 
This will only appear in plots as another point indicating the error of the mean
_intrinsic to the method._ 


```python
# items to be replaced to
replace_with   = ["LOrdered", "fluid", "gel"]
# from those, which have the mark "(1)" at their end
to_replace = [ i+"(1)" for i in replace_with ]
# do the replacement in-place
midata.replace(to_replace=to_replace, value=replace_with, inplace=True)
# show it
midata.loc[midata.particle == "S"].query("`satur index` == 0.0")
```

Good,
now I have a clean data set 
in which I know about some missing values, 
but I also know why and what that it has no effect on what I want to do further.

In case this becomes a problem, I can use _dropna()_ method to skip the values.

description:

```python
# here selecting only particles "T" and grouping by "sterol concentration"
midata.loc[midata.particle=='T'].groupby(by='sterol conc').describe().T
```

```python
midata.loc[:,"starting conf"].value_counts()
```

# Time for transformation!

Nice,
now I have a data set with well defined values in all columns. 

Some columns have ordinal values -
e.g. particle type T - S - R
can be mapped onto 2 - 3 - 4

## turn initial conditions into categorical one-hot enc values
Starting conditions are categorical and
can be assigned using one-hot encoding.

```python
# Get a Pd.Series consisting of all the string categoricals
one_hot_encode_cols = midata.dtypes[midata.dtypes == np.object]  # filtering by string categoricals
# and only those which are type object - eg. particle type of starting cond. are excluded byt the condition in [..]
one_hot_encode_cols = one_hot_encode_cols.index.tolist()  # list of categorical fields

# I wan to encode particles as ordinal values, so let's skip it for now...
one_hot_encode_cols.remove('particle')

midata[one_hot_encode_cols].head(10).T

```

```python
# Do the one hot encoding and create a new dataframe "DF"
df = pd.get_dummies(midata, columns=one_hot_encode_cols)
# J: this is a cool one liner - converts all the conditions into separate columns 
#    containing 0/1 "False/True" vlaues - one hot encoding
df.describe().T
```

Dataframe _df_ now contains one-hot-encoded values of initial condition. Good!

## Optional: R-S-T particles into ordinal values ➡️ now skipped!

Let's encode the particles now:

```python
# this one uses lists/or/tuples,
# which is not very practical for my case.
# documentation:
# https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.OrdinalEncoder.html
ordenc = preprocessing.OrdinalEncoder()
X = [['T', 2], ['S', 3], ['R', 4]]
ordenc.fit(X)
# I keep it here only for completeness.
```

```python
# but then I used also another encoder available, which seems a bit more appropriate:
# Label encoder:
#https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.LabelEncoder.html#sklearn.preprocessing.LabelEncoder
# ple = particle label encoder
ple = preprocessing.LabelEncoder()
# this is not necessary, but I used it during debugging and prototyping:
#categories = ['T', 'S', 'R']
#ple.fit(categories)

#this counts and works auto-magically:
ple.fit(df.particle)
```

applying the transformation to the particle column:

```python
df.particle = ple.transform(df.particle)
df.particle
```

```python
# and this is how I get the original labels back
le.inverse_transform(df.particle)
```

```python
df.info()
```

## Preprocessing finished

Now, I have a data set, which is clean,
contains only numerical values,
for which I have a back-transformation to their original descriptors. 

Let's save it.

```python
df.to_csv("dataFrame_all_sims_processed.csv")
```

# Exploratory plots and additional analysis + processing

This is how the variables perm and compress look in a pair plot (simple visual exploration)

```python
sns.pairplot(df.loc[:, 'compress':'perm'])  #, hue='sterol conc')
```

## Normalization of permeability values
Permeability from simulations does not have a great sense
by itself in the actual values.

In the paper, we use relative values normalized to the permeability at
DOPC. 
That is very similar to the value of POPC, so ...

→ Here, I will use POPC as the reference value for such "scaling". 

```python
# creating a multiindex version of the dataframe
dfmi = df.set_index(keys=df.columns[:3].to_list())
dfmi.head()
```

```python
# group by particle and normalize against permeability saturation index 0.5 == 1
permnorm = dfmi.copy()
norm_label = 'perm'
for particle_ind in ['R', 'S', 'T']:   # use the following if RST converted to ordinal ints: [0, 1, 2]:
    # determine the denominator for scaling:
    norm_denom = float(dfmi[dfmi['starting conf_fluid']==1].loc[pd.IndexSlice[0.5, 0.0, particle_ind], [norm_label]].values)
    print("Norm factor for particle type {} is {}".format(particle_ind, norm_denom))

    # do the transform for the property
    permnorm.loc[pd.IndexSlice[:,:,particle_ind], (norm_label)] = dfmi.loc[pd.IndexSlice[:,:,particle_ind], 
                                                                           (norm_label)].apply(lambda x: x / norm_denom)
    # ... as well as for its error (if exists!)
    try: 
        norm_label_err = norm_label+" err"
        permnorm.loc[pd.IndexSlice[:,:,particle_ind], (norm_label_err)] = dfmi.loc[pd.IndexSlice[:,:,particle_ind], 
                                                                                      (norm_label_err)].apply(lambda x: x / norm_denom)
    except:
        print("Couldn't perform the transform on the adjacent error to the property {}".format(norm_label))
```

Now the comparison between particles is hampered due to the 
additional scaling of their values.
However, that is not that problematic, 
as the differences for liquid-fluid membranes are not large
and they (likely) arise mostly from the slightly different solubility
of the different beads (which is not precisely the same).
That is why there is lower permeability for smaller particles (against expectations). 

This can also be understood as a benefit - 
now the "solubility" is _"scaled away"_ from the trends
and only the size of the permeating particle drives the major changes. 

👉 i.e. we can now simply compare the effects of size 
as the initial permeability (mostly solubility) 
is basically neglected after this normalization. 

```python
# scaling of the permeability values has happened:
permnorm.loc[0.5, 0]
```

Assign the normalized permeabilities back to the *df* Dataframe 
so that the rest of this notebook works without changes:

```python
# assign the normalized permeabilities back to the df Dataframe 
# so that the rest of this notebook works without changes
df = permnorm.reset_index()
```

🚨 Now, the Dataframe **df** contains permeabilities 
**normalized** to the value through POPC separately for each particle type. 
(see discussion above why this is 👌)

## Skew-ness 
Some properties - especially  permeability and compressibility appear very skewed.

That is not surprising, as we know from the measurements, that there are 
order of magnitude differences in the permeability values depending on the 
saturation index and other model parameters. 

The simple method to make the distributions less skewed 
-- and more meaningful revealing the improtant order-of-magnitude changes --
is to use a _log_ transformation.

It looks that both _compress_ as well as the known _perme_ will benefit from the transform let's see...

Here I calculate the skewness, just for the record:

```python
df.skew()
```

```python
# Let's look at what happens to one of these features, when we apply np.log visually.

# which columns appear skewed and will be investigated:
skew_cols = ['perm', 'perm err', 'compress']

# Choose the transform function for the skewed data
trans_func = np.log10

# Choose a field
for field in skew_cols:

    # Create two "subplots" and a "figure" using matplotlib
    fig, (ax_before, ax_after) = plt.subplots(1, 2, figsize=(10, 5))

    # Create a histogram on the "ax_before" subplot
    df[field].hist(ax=ax_before)

    # Apply a log transformation (numpy syntax) to this column
    # does NOT happen in-place
    df[field].apply(trans_func).hist(ax=ax_after)

    # Formatting of titles etc. for each subplot
    ax_before.set(title='before transfrom', ylabel='frequency', xlabel='value')
    ax_after.set(title='after transform', ylabel='frequency', xlabel='value')
    fig.suptitle('Field "{}"'.format(field));
    
    # and print the skew values before/after the same transform
    print("Skew of {} before transform: {}".format(field, df[field].skew()))
    print("Skew of {} after  transform: {}".format(field, df[field].apply(trans_func).skew()))
```

Skeweness of both properties turns to negative values after transformation.

The value for _compress_ is closer to zero, while the value for _perm_ is more skewed than before -
now only negatively. 

Let's see how the pairplot looks after the transform:

```python
sns.pairplot(df.loc[:, ['compress', 'perm', 'perm err']].apply(trans_func))
```

### Perform the log10 transform
.. and save it into a separate dataframe *dft* :

```python
# Perform the skew transformation on a copy of the dataframe:
dft = df.copy()

for col in skew_cols:  
    dft[col] = df[col].apply(trans_func)
    
dft.rename(columns={ col: "log10_"+col for col in skew_cols}, inplace=True)
```

```python
dft.info()
```

## Pair-plot for all particles

Particles are color-coded in the following plot.

If Ordinal encoding (via LabelEncoder, see above) is used for their names,
this is the code:
- R: 0
- S: 1
- T: 2

👇

```python
hue_col = 'particle'
sns.pairplot(dft.loc[:, ['satur index', hue_col, 'sterol conc', 'log10_perm']]  , hue=hue_col)
```

## Pair-plot for particle T 
the tiny-est particle "T" 
yields the most robust resutls 
w.r.t to determinig permeability as 
its size allows for its passing through the membrane
even in the _gel phase_. 

Here I plot the pair-correlations of the properties from various simulations
only for the **T** particle. 

👇

```python
# pair-plot the result:
hue_col = 'sterol conc'
sns.pairplot(dft[dft.particle == 'T'].loc[:, ['satur index', hue_col, 'log10_compress', 'log10_perm']]  , hue=hue_col)
```

## Quick summary of the results - mostly for particle T
The above plots contain valuable information on the effect of sterols and saturation on permeability.
Here summarized as a few bullet points:

- Sterols prevent the membrane from being at the _gel phase_
    - this increases the permeability
- 10% cholesterol can be at a _gel_ phase and has according properties (low permeability)
    - compared to no-cholesterol-DPPC-only, 
      adding just 10% chol. increases permeability mildly (if phase is not affected)
- Sterols gradually decrease permeability from _fluid_ to _LOrdered_
    - this effect is within one order of magnitude
- _gel_ -vs- _fluid_ yield ~2 order of magnitude difference in permeability
- Permeability of a Mixtrue of _fluid_ lipids with 45% sterols 
  has _almost_ one order of magnitude difference from
  Perme. of a mix of DPPC+45% chol.
      - this is in agreement with experimental measurements for water
      - experiments for Fromic acid have some additional chemistry happening 
        (effect of size? - see particles [0,1]→[R,S]


```python
# Let's save the final dataframe including transformation
dft.to_csv("dataFrame_all_sims_processed_log-transf_norm.csv")
```

## Create a dataframe with multiindex for simple cutting and plotting

as the columns in the dataframe are

```python
dft.columns
```

I can simply use the following to make an index from the first columns:

```python
dftmi = dft.set_index(keys=dft.columns[:3].to_list())
dftmi.head()
```

```python
dftmi.loc[0.0, :, 'T'].sort_index()
```

## Plotting dependence of Permeability ⨉ Sterol conc

- for each particle type (color-hue-coded)
- distinguishing different phases POPC/DPPC (top/bottom)
- omitting clearly non-converged results (starting from the incorrect phase and did not make it into the correct)


```python
for satind, init_conf_gel in zip([0.0, 0.5], [True, False]):
    sns.lineplot(x='sterol conc', y='log10_perm', data=dft.query(
        '(`starting conf_gel` == {} or `starting conf_LOrdered` == 1) and `satur index` == {}'.format(init_conf_gel, satind)), 
                 hue='particle').set_title(r'$P_{MD}$ | top: $L_d$ - bottom: $L_o$ resp. $L_\beta$')
    
plt.savefig('plot_perm-log10_sterol-conc_allParticles.png', dpi=150)
```

## Plotting dependence of Permeability ⨉ Satur index

- for each particle type

All starting conditions (and their phases) are plotted in the graph.
They can be easily distinguished by their clear separation.

Obviously, there are convergence problems with the largest Regular particle R (resp. 0).
In the fluid phase Ld, however, all particles 
reveal highly similar behaviour.

```python
for sterolconc in [0, 15, 30, 45]:
    for init_conf_gel in [False, True]:
        sns.lineplot(x='satur index', y='log10_perm', data=dft.query(
            '(`starting conf_gel` == {} or `starting conf_LOrdered` == 1) and `sterol conc` == {}'.
            format(init_conf_gel, sterolconc)), 
            hue='particle').set_title("Sterol concentration: {}%".format(sterolconc))
    plt.savefig('plot_perm-log10_satur-index_allParticles_sterol-conc-{}p.png'.format(sterolconc), dpi=150)
    plt.show()    
```

```python

```
