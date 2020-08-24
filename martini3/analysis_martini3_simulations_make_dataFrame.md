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
import gromacs
import pandas as pd
import os, fnmatch
import matplotlib.pyplot as plt
from IPython.display import display
import MDAnalysis as mda
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
        ir = -1  #till the end of the original string
    return s[il:ir]
```

# Simulation class definition
This class includes all methods and attributes to obtain results from ethanol-membrane intreraction simulations with Martini within MeMBrane project.

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
    column_labels = ("temperature",
                    "sterol type", "sterol conc",
                    "other phosph", "tails", "satur index",
                    "PC conc", "ethanol conc", 
                    "APL", "thickness",
                    "bending", "tilt",
                    "compress")

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
        for fname in [self.top_fname, self.lipidator_fname, self.solvent_density_fname, self.boxx_fname]:
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
            
    
    def calc_permeability(self):
        """calculate permeability from AWH simulation profiles and friction"""
        ============ TBD!
        awh_data  = np.loadtxt(newest_file, comments=("#", "@"))
        fric_data = np.loadtxt(newest_file.replace("awh_", "friction_"), comments=("#", "@"))
        x    = awh_data[:,0]
        fep  = awh_data[:,1]
        fric = fric_data[:,1]
        # OVERRIDING friction data with a constant number
        #fric = np.ones(fep.shape)*2900.0
        # obtain the estimate of deltaG
        dg1   = np.max(fep-fep[-m:].mean())
        dg2   = np.max(fep-fep[:m].mean())
        dgs.append(np.mean([dg1, dg2]))
        dgstds.append(np.std([dg1, dg2]))
        if use_symmetry:
            fepsymm  = symmetrize(fep)
            fricsymm = symmetrize(fric)
            xsymm    = x[:fepsymm[0].shape[0]]
            rxsymm    = np.exp(fepsymm[0]-fepsymm[0][:m].mean())*fricsymm[0]/1000.0    # r(x)
            #rxstdsymm = ( np.exp(fepsymm[0])*fricsymm[1] + np.exp(fepsymm[1])*fricsymm[0] )/1000.0    # error of r(x) through chain rule
            # OVERRIDING error estimate of friction here
            rxstdsymm = ( np.exp(fepsymm[0])*100.0 + np.exp(fepsymm[1])*fricsymm[0] )/1000.0    # error of r(x) through chain rule
            r         = np.trapz(x=xsymm[m:], y=rxsymm[m:]) *2.0   # *2, because i have only half of the profile after symmetrizing it
            rstd      = np.trapz(x=xsymm[m:], y=rxstdsymm[m:]) *2.0
            perm      = 100000.0/r   # in nm/us*e-4 resp cm/s*e-3
            permstd   = perm**2 * rstd /100000.0   # is the same as: perm*rstd/r
            simulations[newest_dirname] = r
            means.append(r)
            stds.append(rstd)
            perms.append(perm)
            permstds.append(permstd)
        else:
            # set the zero-level (solvent) to 0 from one resp. other side
            rx1  = np.exp(fep-fep[-m:].mean())*fric/1000.0    # r(x) ...
            rx2  = np.exp(fep-fep[:m].mean()) *fric/1000.0    # resistance in ns/nm resp s/m
            r = np.trapz(x=x[m:-m], y=np.vstack((rx1, rx2))[:,m:-m])
            simulations[newest_dirname] = r.mean()
            means.append(r.mean())
            stds.append(r.std())
         
    
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
        and compositions (which lipids, concentration of ethanol, what sterols and how much..)
        as instance attributes
        
        This method is HIGHLY SPECIFIC to the ethanol project 
        and the dirname naming conventions therein
        """
        s = self.dirname # just a short-hand abbreviation for the long self.dirname 
        self.sterol_type     = cropstring(s, l="./", r="sterol")
        self.sterol_conc  =int(cropstring(s, l="rol/"+self.sterol_type+"sterol", r="p_sims_data"))
        self.ethanol_conc =int(cropstring(s, "memb_", "pEthanol"))
        self.tails           = cropstring(s, "_tails_", "_heads_")
        self.temperature  =int(cropstring(s, "_Temp_", "K"))
        self.composition_str = cropstring(s, "/0", "_martini")
        # now the phospholipid composition
        self.other_phospholipid = cropstring(s, "_heads_", "-PC_Temp")
        if self.composition_str.startswith("1_"):
            self.pc_conc = 100
        elif self.composition_str.startswith("7_"):
            self.pc_conc = 0
        else:
            r_pc    = int(self.composition_str[2])
            r_other = int(self.composition_str[8])
            self.pc_conc = round(100.0*r_pc/(r_pc+r_other))
        # and finally, the dirname path was parsed, so ...
        self.dirname_parsed  = True
        
    
    @property
    def satur_index(self):
        """
        Assign a saturation index d depending on the used tails
        """
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
            self.temperature
        except:
            self.parse_dirname()
        finally:    
            record = (self.temperature, 
                      self.sterol_type, self.sterol_conc, 
                      self.other_phospholipid, self.tails, self.satur_index, 
                      self.pc_conc, self.ethanol_conc,
                      self.apl, self.thick,
                      self.bend, self.tilt, 
                      self.compress)

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
for f in tpr_files:
    sim = Simulation(os.path.dirname(f))
    records.append(sim.as_record)
```

# Table of simulations

```python
# create large DataFrame data set from records
data = pd.DataFrame.from_records(records, columns=Simulation.column_labels)

# create multiindex from columns 1-8
midata = data.set_index(list(Simulation.column_labels[:8]))
```

```python
# This is how to slice it (an example)
idx = pd.IndexSlice
#          parts of the multi-index ...   , membrane property(ies)
midata.loc[idx[298, "no", 0, "PE", :, 0.0], :].sort_index(level="PC conc")
```

## Test plot and slicing
Here I perform a simple test of plotting and slicing and manipulating of the dataset

```python
pc_concs = midata.index.levels[6]
eth_concs = midata.index.levels[7]

memb_prop = 'tilt'
sterol_type = "no"
sterol_conc = 0
for pltkind in ["kde", ]:
    midata.loc[idx[298, sterol_type, sterol_conc, "PE", :, 0.0], memb_prop].sort_index(
        level="PC conc").plot(kind=pltkind)
    #plt.xlabel("area per lipid [nm]")
    plt.show()
    

for pltkind in ["line", ]:
    for pcc in pc_concs:
        midata.loc[idx[298, sterol_type, sterol_conc, "PE", :, 0.0, pcc, :], memb_prop].sort_index(
            level="PC conc").plot(kind=pltkind, label=pcc)
        plt.legend()
        plt.xticks(range(len(eth_concs)), eth_concs)

```

## Save it!

```python
# CSV is not much larger than other data formats and can be read-in by most software suites including e.g. Excel
midata.to_csv(path_or_buf="dataFrame_all_sims.csv")
```

# Next ...
plotting and working with the dataset will be done in a separate Jupyter notebook, 
so that this one is clean and finishes with just the creation of the dataset.