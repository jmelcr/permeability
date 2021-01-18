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

# Function definitions
a few handy functions with a quite general meaning 

```python
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result
```

```python
def locate(pattern, root=os.curdir):
    '''Locate all files matching supplied filename pattern in and below
    supplied root directory.'''
    for path, dirs, files in os.walk(os.path.abspath(root)):
        for filename in fnmatch.filter(files, pattern):
            yield os.path.join(path, filename)
```

```python
def symmetrize(a):
    """ 
    a - array to be symmetrized mirror-like, with the mirror symmetry placed in the centre
    returns: a tuple with the symmetrized array (half the size), 
             and the resulting standard deviation in the second place of the tuple
    """
    d = a.shape[0]   #assuming index=0
    m = d//2         # the middle-of-array index -- integer division here
    l = a[:m]        # left part
    r = np.ndarray(l.shape)
    for i in range(m):
        r[i] = a[-i-1]
    s = np.average(np.vstack((l, r)), axis=0)
    e = np.std(np.vstack((l, r)), axis=0)
    return (s, e)
```

```python
def cropstring(s, l=None, r=None):
    """
    Returns a cropped substring of a string s
    
     l, r : strings
            strings or patterns (RegEx not implemented!) to match in string s
            for cropping after l and before r
    """
    if isinstance(l, str):
        il = s.find(l) + len(l)
    else:
        il = 0
    if isinstance(r, str):
        ir = s.find(r)
    else:
        ir = len(s)+1  #till the end of the original string
    return s[il:ir]
```

# Simulation class definition
This class includes all methods and attributes to obtain results from ethanol-membrane intreraction simulations with Martini within **permeability** project with Jacopo Fralicciardi and Bert Poolman.

This is prepared in a way that the readout shall be easily concatenated in a Pandas Dataframe for further analysis including multi-indexing etc...

```python
class Simulation:
    '''Simulation class stores and accesses simulation data 
    relevant for ethanol - membrane interaction project with Martini'''
    default_top_fname             = "topol.tpr"
    default_traj_fname            = "traj_comp_pbc.xtc"
    default_lipidator_fname       = "lipidator.out"
    default_solvent_density_fname = "density_solvent.xvg"
    default_boxx_fname            = "boxx.xvg"
    use_symmetry = True
    column_labels = ("satur index",
                     "sterol conc",
                     "particle",
                     "starting conf",
                     "compress",
                     "perm", "perm err")

    def __init__(self, dirname):
        self.dirname = dirname
        # the dirname containing conditions of the simulation
        # is not parsed at the beginning (only done when needed in method parse_dirname)
        self.dirname_parsed = False
        # check presence of required files
        self.top_fname             = Simulation.default_top_fname
        self.traj_fname            = Simulation.default_traj_fname
        self.lipidator_fname       = Simulation.default_lipidator_fname
        self.solvent_density_fname = Simulation.default_solvent_density_fname
        self.boxx_fname            = Simulation.default_boxx_fname
        for fname in [self.top_fname, ]: # self.lipidator_fname, self.solvent_density_fname, self.boxx_fname]:
            try:
                for f in locate(fname, self.dirname):
                    if fname in f:
                        temp_found_file0 = f
                    break
            except:
                print("something went wrong with locating {} in dir:\n{}".format(fname, self.dirname))    
        # get the topology using MDAnalysis
        # this is not used at this moment, so creating mda.Universe is skipped
        if False:
            try:
                self.topology = mda.Universe(os.path.join(self.dirname, self.top_fname)) 
            except:
                raise UserWarning("Error loading topology and trajectory in directory {}".format(self.dirname))
            
            
    
    def calc_perm(self, m=20, use_symmetry=True, const_fric=None):
        """
        calculate permeability from AWH simulation deltaG and friction profiles
        
        parameter m defines the number of samples from each side of the profile
        to use for estimation of the zero-level 
        by taking mean on the left & right (which are obtained independently and have to agree in a converged sim.) 
        and comparing the values
        
        returns tuple with two values (permeability and , its err.estim.)
        """
        
        # locate the awh file (friction file has the same naming convention format, se we will only simply substitute)
        awh_files = locate("awh_*.xvg", self.dirname)
        # choose the newest of the files
        try:
            newest_file = max(awh_files, key=os.path.getctime)
        except:
            print("cannot find newest file in dir {} \n setting perm to None values".format(sim.dirname))
            return (None, None)
        
        # assign the AWH data sets to numpy arrays
        awh_data  = np.loadtxt(newest_file, comments=("#", "@"))
        x    = awh_data[:,0]
        fep  = awh_data[:,1]

        # prepare friction data - from file (implicit) or constant (optional)
        if not const_fric:
            fric_data = np.loadtxt(newest_file.replace("/awh_", "/friction_"), comments=("#", "@"))
            fric = fric_data[:,1]
        else:
            # OVERRIDING friction data with a constant number
            fric = np.ones(fep.shape)*const_fric
            
        # store the AWH, x-coord & friction profile as attributes of this class Simulation
        self.awh_x = x
        self.awh   = fep
        self.fric  = fric
            
        # obtain the estimate of deltaG from left & right
        dg1 = np.max(fep-fep[-m:].mean())
        dg2 = np.max(fep-fep[:m].mean())
        self.dg    = np.mean([dg1, dg2])
        self.dgstd = np.std([dg1, dg2])
        
        # calculate the permeability (through 1/resistance)
        if use_symmetry:
            # apply symmetry rules to the profile and obtain an estimate of error
            fepsymm  = symmetrize(fep)
            fricsymm = symmetrize(fric)
            xsymm    = x[:fepsymm[0].shape[0]]
            rxsymm    = np.exp(fepsymm[0]-fepsymm[0][:m].mean())*fricsymm[0]/1000.0    # r(x) = resistance
            # use value of 100.0 instead of "*fricsymm[1]" to simplify the equation and error estimate?
            rxstdsymm = ( np.exp(fepsymm[0])*fricsymm[1] + np.exp(fepsymm[1])*fricsymm[0] )/1000.0    # error of r(x) through chain rule
            r         = np.trapz(x=xsymm[m:], y=rxsymm[m:]) *2.0   # *2, because i have only half of the profile after symmetrizing it
            rstd      = np.trapz(x=xsymm[m:], y=rxstdsymm[m:]) *2.0
        else:
            # set the zero-level (solvent) to 0 from one resp. other side
            rx1  = np.exp(fep-fep[-m:].mean())*fric/1000.0    # r(x) ...
            rx2  = np.exp(fep-fep[:m].mean()) *fric/1000.0    # resistance in ns/nm resp s/m
            r = np.trapz(x=x[m:-m], y=np.vstack((rx1, rx2))[:,m:-m])

        unit_scaling_factor = 100000.0  # provides unit cm/s *e-3 from simulation units (kT and ps)d
        perm      = unit_scaling_factor/r   # in nm/us*e-4 resp cm/s*e-3
        try: 
            permstd = perm**2 * rstd /unit_scaling_factor   # is the same as: perm*rstd/r
        except:
            permstd = None

        self.__perm_perr = (perm, permstd)
        return (perm, permstd)

    
    @property
    def perm(self):
        try:
            return self.__perm_perr[0]
        except:
            # shall be raised e.g. in case self.apl is not defined yet
            self.calc_perm()  # defines the above self."attribute"
        finally:
            return self.__perm_perr[0]

    @property
    def permerr(self):
        try:
            return self.__perm_perr[1]
        except:
            # shall be raised e.g. in case self.apl is not defined yet
            self.calc_perm()  # defines the above self."attribute"
        finally:
            return self.__perm_perr[1]
        
    
    def read_lipidator(self):
        """get several properties,
        bending and tilt moduli [in units kT and kT/nm^2],
        area per lipid (incl. stderr & stdev) [in units nm^2], 
        from lipidator file"""
        try:
            with open(os.path.join(self.dirname, self.lipidator_fname), "r") as lefile:
                for line in lefile.readlines():
                    linesplit = line.split()
                    if line.startswith("Monolayer Bending"):
                        bending = float(linesplit[linesplit.index('kappa')+1])
                    if line.startswith("Bilayer Tilt"):
                        tilt    = float(linesplit[linesplit.index('kappa')+1])
                    if line.startswith(" Area AVG"):
                        apl       = float(linesplit[linesplit.index('[bilayer]:')+1]) *0.01
                        aplstderr = float(linesplit[linesplit.index('[bilayer]:')+4]) *0.01
                        aplstdev  = float(linesplit[linesplit.index('[bilayer]:')+6]) *0.01
                        
        except IOError:
            print("Problems loading the lipidator file in dir {}".format(self.dirname))
        except:
            print("Generally some problems in loading Lipidator output -- check code.")
        finally:
            try:
                self.__bend       = bending
                self.__tilt       = tilt
                self.__apl        = apl
                self.__apl_stderr = aplstderr
                self.__apl_stdev  = aplstdev
            except:
                # will be raised e.g. in case any of the above is not defined 
                # (e.g. weren't accessed via the method above)
                self.__bend = self.__tilt = self.__apl = self.__apl_stderr = self.__apl_stdev = None
        return (bending, tilt, apl, aplstderr, aplstdev)
    
    
    @property
    def apl(self):
        try:
            return self.__apl
        except:
            # shall be raised e.g. in case self.apl is not defined yet
            self.read_lipidator()  # defines the above self.apl
        finally:
            return self.__apl
    

    @property
    def tilt(self):
        try:
            return self.__tilt
        except:
            # shall be raised e.g. in case self.apl is not defined yet
            self.read_lipidator()  # defines the above self.apl
        finally:
            return self.__tilt
    

    @property
    def bend(self):
        try:
            return self.__bend
        except:
            # shall be raised e.g. in case self.apl is not defined yet
            self.read_lipidator()  # defines the above self.apl
        finally:
            return self.__bend
    
    
    def calc_thick(self):
        """
        method calculating Luzzati thickness (hydrophobic thickness) 
        from a density profile xvg file
        """
        try:
            # variable wat also accounts for generally any solvent (e.g. water+ethanol solution)
            dens_arr = np.loadtxt(os.path.join(self.dirname, self.solvent_density_fname), comments=["#", "@"])
            wat = dens_arr[:,1]
            x   = dens_arr[:,0]
            d_z = np.abs(x[-1]-x[0])         # simulation box z-length , aka repeat distance
            bulkwatdens = np.mean(wat[:10])
            wat /= bulkwatdens
            h   = d_z - np.trapz(wat, x)
            #print(x, wat, bulkwatdens, d_z, h)
        except IOError:
            print("Problems loading the density_ethanol_water XVG file in dir {}".format(self.dirname))
            h = None
        except:
            print("Generally some problems in the density_ethanol_water loading resp. Luzzati thickness calc -- check code.")
            h = None
        finally:
            self.__thick = h
            return thick                


    @property
    def thick(self):
        try:
            return self.__thick
        except:
            # shall be raised e.g. in case self.__thick is not defined yet
            self.calc_thick()
        finally:
            return self.__thick

            
    #    Copied from below - turn into methods!
    def calc_compressibility(self):
        """
        method calculating 
        Area Compressibility modulus [in units kT/nm^2]
        from Box-x xvg file - a time-series of box-x (data : np.array)
        using the fluctuation theorem  
        """
        try:
            data_arr = np.loadtxt(os.path.join(self.dirname, self.boxx_fname), comments=["#", "@"])
            x_t = data_arr[:,1]
            a_t = x_t**2
            a_0 = a_t.mean()
            var =  a_t.var()
            ka = a_0/var
        except IOError:
            print("Problems loading the box-x XVG file in dir {}".format(self.dirname))
            ka = None
        except:
            print("Generally some problems in the Box-x loading/calc -- check code.")
            ka = None
        finally:
            self.__compress = ka
            return ka
           
            
    @property
    def compress(self):
        try:
            return self.__compress
        except:
            # shall be raised e.g. in case self.apl is not defined yet
            self.calc_compressibility()  # defines the above self."attribute"
        finally:
            return self.__compress
 

    def parse_dirname(self):
        """
        Dirname parser assigning values to 
        simulation conditions (temperature)
        and compositions (saturation number d, which lipids, what sterols and how much..)
        as instance attributes
        
        This method is HIGHLY SPECIFIC to the permeability project 
        and the dirname naming conventions therein
        """
        s = self.dirname # just a short-hand abbreviation for the long self.dirname 
        self.sterol_type     = cropstring(s, l="./", r="sterol")
        self.sterol_conc  =int(cropstring(s, l=self.sterol_type+"sterol", r="p_PC"))
        self.temperature  = 293  # as in the main mdp file
        self.particle = cropstring(s, "perm_particle_EOL", "/sat-unsat")
        self.d = float(cropstring(s, "_d_", "/sim1"))/100.0  # saturation index d
        self.starting_conf = cropstring(s, "_awh_")
        # and finally, the dirname path was parsed, so ...
        self.dirname_parsed  = True
        
    
    @property
    def satur_index(self):
        """
        Assign a saturation index d depending on the used tails
        
        NOT USED in this project
        """
        raise Warning("Invoking function >> satur_index <<, which shall not be used!")
        try:
            self.tails
        except:
            self.parse_dirname()
        finally:
            t = self.tails
            if t == "DP":
                d = 0.0
            elif t == "PO":
                d = 0.5
            elif t == "DO" or t == "PI":
                d = 1.0
            elif t == "DI":
                d = 2.0
            else:
                d = None
        self.sat_ind = d
        return d
    
    
    @property
    def as_record(self):
        """
        Export the single Simulation instance as a record (a tuple of items)
        The order is given by the Class variable Simulation.column_labels
        """
        try:
            self.d
        except:
            self.parse_dirname()
        finally:    
            record = (self.d, 
                      self.sterol_conc, 
                      self.particle,
                      self.starting_conf,
                      self.compress,
                      self.perm, self.permerr)

        return record
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
midata.loc[idx[:, 0, "S"], :].sort_index(level="satur index").head()
```

## Test plot and slicing
Here I perform a simple test of plotting and slicing and manipulating of the dataset

```python
sns.pairplot(data.loc[:, ['particle', 'perm', 'satur index', 'sterol conc']], hue="particle")
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
    if "EOLS" in s.dirname and not "(1)" in s.dirname:
        print(s.dirname)
        try:
            s.parse_dirname()
            print(s.d, s.sterol_conc, s.particle)
        except:
            pass
        if s.d == 0.5:
            po_sims.append(s)
            print(" ➡ added this one to the 💧 PO-list\n")
        if s.d == 0.0:
            dp_sims.append(s)
            print(" ➡ added this one to the 🧊 DP-list\n")
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
def prep_to_plot(x, y, shift_to_zero=True):
    """Function symmetrizes and shifts the profiles for plotting."""
    ysymm_arr = symmetrize(y)
    ys    = ysymm_arr[:][0]
    yserr = ysymm_arr[:][1]
    x_half = -x[:len(ys)]
    if shift_to_zero:
        delta = ys[:50].mean()
        ys -= delta
    return (x_half, ys, yserr)
```

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
        if (xxpc == 'POPC' and s.starting_conf != 'gel') or (xxpc == 'DPPC' and s.starting_conf != 'fluid'):
            (x_half, fricsym, fricsym_err) = prep_to_plot(s.awh_x, s.fric, shift_to_zero=False)
            x_half = -s.awh_x[:len(fricsym)]
            plt.plot(x_half, fricsym)
            plt.fill_between(x=x_half, 
                             y1=fricsym+fricsym_err,
                             y2=fricsym-fricsym_err,
                             label="{}, {}% {}sterol".format(xxpc, s.sterol_conc, s.sterol_type),
                             alpha=0.9)
    plt.legend()
    plt.ylabel("friction / ps nm$^{-2}$ (kT)$^{-1}$")
    plt.xlabel("distance / nm")
    plt.savefig("friction_profiles_{}_sterol-concs-EOLS.png".format(xxpc), dpi=150, bbox_inces='tight')
    plt.show()
```

```python

```
