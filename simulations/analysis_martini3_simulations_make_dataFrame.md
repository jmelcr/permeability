---
jupyter:
  jupytext:
    formats: ipynb,py:percent,md
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
import numpy as np
#import gromacs
import pandas as pd
import os, fnmatch
import matplotlib.pyplot as plt
import seaborn as sns
from IPython.display import display
import MDAnalysis as mda
import pickle
%pylab inline
```

# Simulation class definition
This class includes all methods and attributes to obtain results from ethanol-membrane intreraction simulations with Martini within **permeability** project with Jacopo Fralicciardi and Bert Poolman.

This is prepared in a way that the readout shall be easily concatenated in a Pandas Dataframe for further analysis including multi-indexing etc...


### Function definitions
a few handy functions with a quite general meaning 

```python
# import the definitions from the separate script
from permeana_class_simulation import *
```

# Search through the simulations
on my local disk on md39.

In the following section, I will loop through the simulation folders and 
obtain results from my analysis scripts.

This will be put together in a Pandas Dataframe with multiindexing at the end. 

The multiindex will use properties as parsed from the directory names.

```python
tpr_files = find("topol.tpr", os.curdir)
```

```python
records = []
sims    = []
for f in tpr_files:
    sim = Simulation(os.path.dirname(f))
    sims.append(sim)
    if (not "prep" in sim.dirname) and ("titration" in sim.dirname) :
        try:
            records.append(sim.as_record)
        except:
            print("Troubles loading simulation in dir {} \n Skipping it.".format(sim.dirname))
```

# Table of simulations

```python
# create large DataFrame data set from records
data = pd.DataFrame.from_records(records, columns=Simulation.column_labels)

data.head()
```

```python
# create multiindex from columns 1-8
midata = data.set_index(list(Simulation.column_labels[:4]))

# This is how to slice it (an example)
idx = pd.IndexSlice
#          parts of the multi-index ...   , membrane property(ies)
midata.loc[idx[:, 0, "S"], :].sort_index(level="satur index").head(6)
```

## Test plot and slicing
Here I perform a simple test of plotting and slicing and manipulating of the dataset

```python
sns.pairplot(data.loc[:, ['particle', 'perm', 'satur index', 'sterol conc', 'thickness', 'compress']], hue="particle")
```


## Save it!

```python
# Save the dataframe with the analyzed trajectories
# CSV is not much larger than other data formats and can be read-in by most software suites including e.g. Excel
data.to_csv(path_or_buf="dataFrame_all_sims.csv")

# Save the list of simulation objects using pickle
with open("objects_all_sims.pickle", "wb") as f: pickle.Pickler(f).dump(sims)
```

# Dataframe analysis separately
plotting and working with the Pandas Dataframe dataset will be done in a separate Jupyter notebook, 
so that this one is clean and finishes with just the creation of the dataset.

# AWH profiles plotting
As working with the pickled objects also requires defining the Simulation class,
I will make the AWH-profile plots here in this notebook. 

## PO– & DP– plots
In this plot, I want to show the profiles of 
simulations with 0,15,30,45% cholesterol
. + DOPC as a reference profile

Here, I find and select the corresponding simulations:

```python
po_sims = []
dp_sims = []
for s in sims:
    # using only Tiny particle as permeant
    if "EOLS" in s.dirname and not "(1)" in s.dirname and not "AWHk" in s.dirname:
        print(s.dirname)
        try:
            s.parse_dirname()
            print(s.d, s.sterol_conc, s.particle)
        except:
            pass
        try:
            if s.d == 0.5:
                po_sims.append(s)
                print(" ➡ added this one to the 💧 PO-list\n")
            if s.d == 0.0:
                dp_sims.append(s)
                print(" ➡ added this one to the 🧊 DP-list\n")
        except:
            pass

    if "DOPC" in s.dirname:
        print("➡️ Found DOPC simulation in dir {}".format(s.dirname))
        dopc_sim = s
        
```

Now, sort the lists after sterol concentrations:

```python
for sim_list in [po_sims, dp_sims]:
    sim_list.sort(key=lambda s: s.sterol_conc)
```

### Free energy profiles

in the code below,
I plot the symmetrized AWH-free energy profiles.
The differences from symmetrization 
give rise to the error estimates included in the plot. 

The plot can be futher improved 
by using a simulation snapshot as a background. 

```python
for sim_list, xxpc in zip([po_sims, dp_sims], ["POPC", "DPPC"]):
    # DOPC - reference - at the background
    s = dopc_sim
    try:
        (x_half, awhsym, awhsym_err) = prep_to_plot(s.awh_x, s.awh)
        plt.plot(x_half, awhsym, color='black')
        plt.fill_between(x=x_half, 
                         y1=awhsym+awhsym_err,
                         y2=awhsym-awhsym_err,
                         label="DOPC, {}% {}sterol".format(s.sterol_conc, s.sterol_type),
                         color='black', alpha=0.7)
    except:
        print("troubles with simulation in {}".format(s.dirname))
    
    for s in sim_list:
        if (xxpc == 'POPC' and s.starting_conf != 'gel') or (xxpc == 'DPPC' and s.starting_conf != 'fluid'):
            try:
                    (x_half, awhsym, awhsym_err) = prep_to_plot(s.awh_x, s.awh)
                    plt.plot(x_half, awhsym)
                    plt.fill_between(x=x_half, 
                                     y1=awhsym+awhsym_err,
                                     y2=awhsym-awhsym_err,
                                     label="{}, {}% {}sterol".format(xxpc, s.sterol_conc, s.sterol_type))
            except:
                print("troubles plotting simulation in {}".format(s.dirname))

    plt.legend()
    plt.ylim([-1, 11])
    plt.ylabel("Free energy / kT")
    plt.xlabel("distance / nm")
    plt.savefig("awh_dG_profiles_{}_sterol-concs-EOLS.png".format(xxpc), dpi=150, bbox_inces='tight')
    plt.show()
```

### Friction profiles

in the code below,
I plot the symmetrized friction profiles.
The differences from symmetrization 
give rise to the error estimates included in the plot. 

These plots would benefit from further processing,
namely smoothening of the noise. 

However, at the current status it makes the necessary points. 

```python
for sim_list, xxpc in zip([po_sims, dp_sims], ["POPC", "DPPC"]):
    for s in sim_list:
        try:
            if (xxpc == 'POPC' and s.starting_conf != 'gel') or (xxpc == 'DPPC' and s.starting_conf != 'fluid'):
                (x_half, fricsym, fricsym_err) = prep_to_plot(s.awh_x, s.fric, shift_to_zero=False, filt_freq=0.150)
                x_half = -s.awh_x[:len(fricsym)]
                plt.plot(x_half, fricsym)
                plt.fill_between(x=x_half, 
                                 y1=fricsym+fricsym_err,
                                 y2=fricsym-fricsym_err,
                                 label="{}, {}% {}sterol".format(xxpc, s.sterol_conc, s.sterol_type),
                                 alpha=0.9)
        except:
            pass
    plt.legend()
    plt.ylabel("friction / ps nm$^{-2}$ (kT)$^{-1}$")
    plt.xlabel("distance / nm")
    plt.savefig("friction_profiles_{}_sterol-concs-EOLS.png".format(xxpc), dpi=150, bbox_inces='tight')
    plt.show()
```

```python
for sim_list, xxpc in zip([po_sims, dp_sims], ["POPC", "DPPC"]):
    for s in sim_list:
        try:
            if (xxpc == 'POPC' and s.starting_conf != 'gel') or (xxpc == 'DPPC' and s.starting_conf != 'fluid'):
                (x_half, fricsym, fricsym_err) = prep_to_plot(s.awh_x, s.fric, 
                                                              shift_to_zero=False, 
                                                              filt_freq=0.150)
                x_half = -s.awh_x[:len(fricsym)]
                plt.plot(x_half, fricsym)
                plt.fill_between(x=x_half, 
                                 y1=fricsym+fricsym_err,
                                 y2=fricsym-fricsym_err,
                                 label="{}, {}% {}sterol".format(xxpc, s.sterol_conc, s.sterol_type),
                                 alpha=0.9)
        except:
            pass
    plt.legend()
    plt.yscale("log")
    plt.ylim([9e2, 5e5])
    plt.ylabel("friction / ps nm$^{-2}$ (kT)$^{-1}$")
    plt.xlabel("distance / nm")
    plt.savefig("friction_profiles_{}_sterol-concs-EOLS-log.png".format(xxpc), dpi=150, bbox_inces='tight')
    plt.show()
```

```python

```
