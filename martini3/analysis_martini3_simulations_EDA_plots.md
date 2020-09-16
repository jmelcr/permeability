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
from cycler import cycler
from sklearn.decomposition import PCA
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
midata.describe()
# here selecting only particles "S" and grouping by "sterol concentration"
midata.loc[midata.particle=='S'].groupby(by='sterol conc').describe().T
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

Starting conditions are categorical and
can be assigned using one-hot encoding.

```python
# Get a Pd.Series consisting of all the string categoricals
one_hot_encode_cols = midata.dtypes[midata.dtypes == np.object]  # filtering by string categoricals
# and only those which are type object - eg. particle type of starting cond. are excluded byt the condition in [..]
one_hot_encode_cols = one_hot_encode_cols.index.tolist()  # list of categorical fields

# I wan to encode particles as ordinal values, so let's skip it for now...
one_hot_encode_cols.remove('particle')

midata[one_hot_encode_cols].head().T

```

```python
# Do the one hot encoding and create a new dataframe "DF"
df = pd.get_dummies(midata, columns=one_hot_encode_cols)
# J: this is a cool one liner - converts all the conditions into separate columns 
#    containing 0/1 "False/True" vlaues - one hot encoding
df.describe().T
```

Dataframe _df_ now contains one-hot-encoded values of initial condition. Good!

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

### Preprocessing finished

Now, I have a data set, which is clean,
contains only numerical values,
for which I have a back-transformation to their original descriptors. 

Let's save it.

```python
df.to_csv("dataFrame_all_sims_processed.csv")
```

```python
sns.pairplot(df, hue='sterol conc')
```
