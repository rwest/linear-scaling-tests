
# coding: utf-8

# # Modeling dry reforming of methane using RMG-Cat
# Based on the notebook to accompany the manuscript: <br/>
# *"Automatic generation of microkinetic mechanisms for heterogeneous catalysis"* by<br/>
# C. Franklin Goldsmith, School of Engineering, Brown University, franklin_goldsmith@brown.edu, and<br/>
# Richard H. West, Department of Chemical Engineering, Northeastern University, r.west@northeastern.edu
# 
# To demonstrate the capabilities of automatic mechanism generation for heterogeneous catalysis we apply our mechanism generator software RMG-Cat to the problem of methane dry reforming on nickel. Comparison is made to the mechanism developed by Olaf Deutschmann and coworkers (Delgado, K.; Maier, L.; Tischer, S.; Zellner, A.; Stotz, H.; Deutschmann, O. Catalysts 2015, 5, 871–904.)

# First, we print what git commit we were on when we ran this notebook, for both the source code (RMG-Py) and the database. 

# In[1]:

get_ipython().run_cell_magic(u'bash', u'', u'cd $RMGpy\npwd\ngit log -n1 --pretty=oneline\ncd ../RMG-database\npwd\ngit log -n1 --pretty=oneline')


# ## Model generation
# We start with a base input file to generate a mechanism for CH4 plus CO2.
# First we print the input file we'll use to generate the model.

# In[2]:

get_ipython().magic(u'cat base/input.py')


# Then we try running it. This will take a couple of minutes.

# In[3]:

get_ipython().run_cell_magic(u'bash', u'', u'python $RMGpy/rmg.py base/input.py > /dev/null\ntail -n12 base/RMG.log')


# There are 52 species and 135 reactions.

# ## Data processing
# Next we will import some libraries and set things up to start importing and analyzing the simulation results.

# In[4]:

get_ipython().magic(u'matplotlib inline')
from matplotlib import pyplot as plt
import matplotlib

# The default output Type 3 (Type3) fonts can't be edited in Adobe Illustrator
# but Type 42 (TrueType) fonts can be, making it easier to move labels around
# slightly to improve layout before publication.
matplotlib.rcParams['pdf.fonttype'] = 42 

# Seaborn helps make matplotlib graphics nicer
import seaborn as sns
sns.set_style("ticks")
sns.set_context("paper", font_scale=1.5, rc={"lines.linewidth": 2.0})

import os
import re
import pandas as pd
import numpy as np
import shutil
import subprocess


# In[5]:

def get_last_csv_file(job_directory):
    """
    Find the CSV file from the largest model in the provided job directory.
    
    For CSV files named `simulation_1_13.csv` you want 13 to be the highest number.
    """
    solver_directory = os.path.join(job_directory,'solver')
    csv_files = sorted([f for f in os.listdir(solver_directory) if f.endswith('.csv') ],
                       key=lambda f: int(f[:-4].split('_')[2]))
    return os.path.join(solver_directory, csv_files[-1])
    
job_directory = 'base'
get_last_csv_file(job_directory)


# We will use Pandas to import the csv file

# In[6]:

last_csv_file = get_last_csv_file(job_directory)
data = pd.read_csv(last_csv_file)
data


# In[7]:

def get_pandas_data(job_directory):
    """
    Get the last CSV file from the provided job directory,
    import it into a Pandas data frame, and tidy it up a bit.
    """
    last_csv_file = get_last_csv_file(job_directory)
    data = pd.read_csv(last_csv_file)
    
    # Make the Time into an index and remove it as a column
    data.index = data['Time (s)']
    del data['Time (s)']
    # Remove inerts that RMG added automatically but we're not using
    for i in 'Ar He Ne'.split():
        del data[i]
    # Remove the Volume column
    del data['Volume (m^3)']
    # Set any amounts that are slightly negative (within the ODE solver's ATOL tolerance), equal to zero
    # to allow 'area' plots without warnings or errors.
    # Anything more negative than -1e-16 probably indicates a bug in RMG and should not be hidden here ;-)
    data[(data<0) & (data>-1e-16)] = 0
    return data


# In[8]:

def rename_columns(data):
    """
    Removes the number (##) from the end of the column names, in place,
    unless doing so would make the names collide.
    Also renames a few species so the plot labels match the names in the manuscript.
    """
    import re
    old = data.columns
    new = [re.sub('\(\d+\)$','',n) for n in old]
    # don't translate names that would no longer be unique
    mapping = {k:v for k,v in zip(old,new) if new.count(v)==1}
    data.rename(columns=mapping, inplace=True)
    
    # Now change a few species that are named differently in the manuscript
    # than in the thermodynamics database used to build the model,
    # so that the plot labels match the manuscript.
    mapping = {'COX':'COvdwX', 'OCX': 'CO=X', 'C2H3X':'CH3CX', 'C2H3OX':'CH3CXO'}
    data.rename(columns=mapping, inplace=True)


# In[9]:

data1 = get_pandas_data('base')
rename_columns(data1)
data1.columns


# In[10]:

# Test it with some plots
data1[['CH4', 'CO2']].plot.line()
data1[['CO', 'H2']].plot.line()
data1[['H2O']].plot.line()


# In[11]:

species_names = data1.columns
# just the gas phase species that aren't always zero:
gas_phase = [n for n in species_names if 'X' not in n and (data1[n]>0).any()]
# all the surface species
surface_phase = [n for n in species_names if 'X' in n]
surface_phase.remove('X')
data1[gas_phase].plot.line()
data1[surface_phase].plot.line()


# In[12]:

print "Significant species (those that exceed 0.001 mol at some point)"
significant = [n for n in data1.columns if(data1[n]>0.001).any()]
with sns.color_palette("hls", len(significant)):
    data1[significant].plot.area(legend='reverse')


# In[13]:

surface = [n for n in data1.columns if 'X' in n and n!='X' and (data1[n]>1e-6).any() ]
print "The {} surface species that exceed 1e-6 mol at some point".format(len(surface))
total_sites = max(data1['X'])
with sns.color_palette('Set3',len(surface)):
    (data1[surface]/total_sites).plot.area(legend='reverse')
    plt.ylabel('site fraction')


# # Effect of binding energies

# In[124]:

# First, make a series of input files in separate directories

with open(os.path.join('base', 'input.py')) as infile:
    input_file = infile.read()
    
base_directory = 'binding_energies'
def directory(carbon,oxygen):
    return os.path.join(base_directory, "c{:.3f}o{:.3f}".format(carbon,oxygen))

def make_input(binding_energies):
    """
    Make an input file for the given (carbon,oxygen) tuple (or iterable) of binding energies
    and return the name of the directory in which it is saved.
    """
    carbon, oxygen = binding_energies
    output = input_file
    out_dir = directory(carbon, oxygen)
    carbon_string = "'C':({:f}, 'eV/molecule')".format(carbon)
    output = re.sub("'C':\(.*?, 'eV/molecule'\)", carbon_string, output)
    oxygen_string = "'O':({:f}, 'eV/molecule')".format(oxygen)
    output = re.sub("'O':\(.*?, 'eV/molecule'\)", oxygen_string, output)
    os.path.exists(out_dir) or os.makedirs(out_dir)
    out_file = os.path.join(out_dir, 'input.py')
    with open(out_file,'w') as outfile:
        outfile.write(output)
    shutil.copy(os.path.join('base','run.sh'), out_dir)
    return out_dir

    
print make_input(-8,-3.5)
    


# In[137]:

def run_simulation(carbon, oxygen):
    """
    Assuming a job file already exists, run it.  This one is local.
    """
    job_directory = directory(carbon, oxygen)
    print job_directory
    assert os.path.exists(job_directory)
    return subprocess.check_call('./run.sh', cwd=job_directory)


# In[138]:

make_input(-8,-3.5)
run_simulation(-8, -3.5)


# In[135]:




# In[92]:

directory(-8,-3.5)


# In[79]:

carbon_range = (-8, -3)
oxygen_range = (-3.5, 0)
grid_size = 5
mesh  = np.mgrid[carbon_range[0]:carbon_range[1]:grid_size*1j, oxygen_range[0]:oxygen_range[1]:grid_size*1j]
mesh


# In[80]:

experiments = mesh.reshape((2,-1)).T
experiments


# In[156]:

carbons, oxygens  = experiments.T
carbons, oxygens


# In[157]:

map(make_input, carbons, oxygens)


# In[158]:

import multiprocessing
pool = multiprocessing.Pool()
pool.map_async(make_input, carbons, oxygens)


# In[81]:

answer = map(sum,experiments)


# In[82]:

rate = np.reshape(answer, (5,5))


# In[83]:

ax = sns.heatmap(rate)


# In[86]:

plt.imshow(rate, interpolation='spline16', origin='lower', extent=(-8,-3, -3.5,0), aspect='equal')


# # STOP HERE.
# stuff below is left over from old notebook

# In[14]:

raise NotImplementedError("Stop here.")


# 
# ## Model generation: with RMG-Cat reaction families
# Now lets look at the version that generates a mechanism by applying reaction families.
# First, inspect how the input file differs from the one above. 

# In[13]:

get_ipython().run_cell_magic(u'bash', u'', u'python $RMGpy/rmg.py ch4_co2_families/input.py > /dev/null\ntail -n12 ch4_co2_families/RMG.log')


# ## Data processing
# First, we make a few plots of the new model.
# Mostly just to show how it can be done, and to see what the results look like. These aren't very pretty as they're not going in the manuscript.

# In[14]:

data2 = get_pandas_data('ch4_co2_families')
rename_columns(data2)
data2.columns


# In[15]:

get_last_csv_file('ch4_co2_families')


# In[16]:

data2[['CH4', 'CO2']].plot.line()
data2[['CO', 'H2']].plot.line()
data2[['H2O']].plot.line()


# In[17]:

print "All species"
data2.plot.area()
plt.show()


# In[18]:

print "Significant species (those that exceed 0.001 mol at some point)"
significant = [n for n in data2.columns if(data2[n]>0.001).any()]
with sns.color_palette("hls", len(significant)):
    data2[significant].plot.area(legend='reverse')


# In[19]:

surface = [n for n in data2.columns if 'X' in n and n!='X' and (data2[n]>1e-6).any()]
print "The {} surface species that exceed 1e-6 mol at some point".format(len(surface))
with sns.color_palette('Set3',len(surface)):
    data2[surface].plot.area(legend='reverse')


# In[20]:

species_names = data2.columns
gas_phase = [n for n in species_names if 'X' not in n and (data2[n]>0).any()]
surface_phase = [n for n in species_names if 'X' in n]
surface_phase.remove('X')
print "Total moles of gas"
data2[gas_phase].sum(axis=1).plot()
plt.show()
print "All gas phase species"
data2[gas_phase].plot.line()
plt.show()
print "All surface species"
data2[surface_phase].plot.line()
plt.show()


# In[21]:

print "A comparison of the two reactants across the two models"
ax = plt.subplot()
data1[['CH4', 'CO2']].plot.line(ax=ax) # seed
data2[['CH4', 'CO2']].plot.line(ax=ax, style=':') # families
plt.show()


# ## Model comparison
# Now we will start making some prettier plots for the manuscript, comparing the two models

# In[22]:

plt.rcParams.update({'mathtext.default':  'regular' }) # make the LaTeX subscripts (eg. CH4) use regular font
sns.set_context("paper", font_scale=1, rc={"lines.linewidth": 1.5}) # Tweak the font size and default line widths

def comparison_plot(subplot_axis=None, species='CH4'):
    label = re.sub('(\d)',r'$_\1$',species)
    ax1 = subplot_axis or plt.subplot()
    
    plt.locator_params(nbins=4) # fewer tick marks
    
    ax1.plot(0, data1[species].iloc[0], 'ko')
    ax1.plot(data1.index, data1[species],'k--', linewidth=2.5)
    ax1.plot(data2.index, data2[species],'r-')
    plt.xlim(0,1)
    ax1.set_ylabel('moles')
    ax1.set_xlabel('time (s)')
    
    # Legend
    dummy = matplotlib.patches.Patch(alpha=0)
    plt.legend([dummy],[label], loc='best', fontsize=12)
    
    plt.tight_layout()

fig = plt.figure(figsize=(2.5,2.5))
ax1 = comparison_plot()


# In[23]:

fig = plt.figure(figsize=(7,4))
for n, species in enumerate('CH4 CO H2O CO2 H2'.split()):
    ax = plt.subplot(2,3,n+1)
    comparison_plot(ax, species)

ax = plt.subplot(2,3,6)
ax.axis('off')
red = matplotlib.lines.Line2D([], [], color='r', label='RMG-Cat')
dash = matplotlib.lines.Line2D([], [], color='k', linestyle='--', linewidth=2.5, label='Delgado et al.')

plt.legend(handles=[red, dash], loc='best',fontsize=12)
plt.tight_layout()
plt.savefig('Multi-panel gas comparison.pdf', bbox_inches='tight')


# In[24]:

def extras_plot(subplot_axis=None, species_list=['C2H6','CH2O'], label_positions=None):
    """
    Plot the requested species on one plot,
    with just the 'data2' (RMG-Cat) values.
    Useful for species not in the Delgado model.
    """
    ax1 = subplot_axis or plt.subplot()
    plt.locator_params(nbins=4) # fewer tick marks
    plt.ticklabel_format(style='sci', axis='y', scilimits=(3,3))
    
    for i,species in enumerate(species_list):
        label = re.sub('(\d)',r'$_\1$',species)
        ax1.semilogy(data2.index, data2[species],'-', label=label)
        plt.ylim(ymin=1e-9)
        plt.xlim(0,1)
        ax1.set_ylabel('moles')
        ax1.set_xlabel('time (s)')
        # Manually constructed label
        x = label_positions[species] if label_positions and species in label_positions else 0.5
        y = data2[x:][species].iloc[0]
        ax1.text(x,y,label,
                 verticalalignment='bottom',
                 color=sns.color_palette()[i],
                 fontsize=10)
        plt.tight_layout()

# test it
fig = plt.figure(figsize=(4,4))
with sns.color_palette("Dark2", 3):
    extras_plot(species_list='C2H6 CH2O C3H8'.split())


# In[25]:

# Everything that exceeds some lower threshold
gas_phase = [n for n in species_names if 'X' not in n and (data2[n]>1.1e-9).any()]
# Remove things already in the comparison plots
gas_phase = [n for n in gas_phase if n not in 'CH4 CO H2O CO2 H2 N2'.split()]
print "Extra gas phase species of note are {}.".format(", ".join(gas_phase))
fig = plt.figure(figsize=(3.25,3.25))
with sns.color_palette("Dark2", len(gas_phase)):
    extras_plot(species_list=gas_phase, 
               label_positions={'C2H5':0.15, 'C3H8':0.1, 'CH3OH':0.7, 'H':0.3})
plt.savefig('Gas extra species (many).pdf',bbox_inches='tight')


# In[26]:

def multi_comparison_plot(subplot_axis=None, species_list=['X', 'HX'], label_positions=None):
    """
    Plot many surface species AND their comparison (if there is one) from the Delgado model
    on a semilog plot. 
    """
    ax1 = subplot_axis or plt.subplot()
    plt.locator_params(nbins=4) # fewer tick marks
    for i,species in enumerate(species_list):
        assert 'X' in species, "For surface species only (y axis is normalized for fractional coverage)."
        label = re.sub('(\d)',r'$_\1$',species)
        label = label.replace('X','*')
        if label=='*': label = 'vacant'
        if species in data1.columns:
            sites = data1['X'].iloc[0]
            ax1.semilogy(data1.index, data1[species]/sites,'k--', linewidth=2.5)
        sites = data2['X'].iloc[0]
        ax1.semilogy(data2.index, data2[species]/sites,'-')
        plt.xlim(0,1)
        plt.ylim(ymin=1e-6)
        ax1.set_ylabel('site fraction')
        ax1.set_xlabel('time (s)')
        ## Manually constructed label
        x = label_positions[species] if label_positions and species in label_positions else 0.5
        y = data2[x:][species].iloc[0] / sites
        ax1.text(x,y,label, 
                 fontsize=10,
                 verticalalignment='bottom',
                 color=sns.color_palette()[i])

# test it
multi_comparison_plot()


# In[27]:

species_names = data2.columns
gas_phase = [n for n in species_names if 'X' not in n and (data2[n]>0).any()]

sites = data2['X'].iloc[0] # initial number of suface sites
# get only adsorbates that exceed a site fraction of 1e-6
surface_phase = [n for n in species_names if 'X' in n and (data2[n]>sites*1e-6).any()]
surface_phase


# In[28]:

label_positions = {'HX':0.4, 'COX':0.6, 'CHX':0.7, 'OX':0.4, 
                   'CHOX':0.6, 'CX':0.3, 'CH2X':0.05,
                   'C2HOX':0.8, 'CHO2X':0.7}
fig = plt.figure(figsize=(5,5))
with sns.color_palette("Dark2", len(surface_phase)):
    ax1 = multi_comparison_plot(species_list=surface_phase, label_positions=label_positions)
plt.tight_layout()
plt.savefig('Surface comparison semilog.pdf', bbox_inches='tight')


# In[29]:

def multi_comparison_loglogplot(subplot_axis=None, species_list=['X', 'HX'], label_positions=None):
    """
    Plot many things AND their comparison (if there is one) from the Delgado model
    on a log-log plot
    """
    ax1 = subplot_axis or plt.subplot()
    for i,species in enumerate(species_list):
        label = re.sub('(\d)',r'$_\1$',species)
        label = label.replace('X','*')
        if label=='*': label='vacant'
        if species in data1.columns:
            sites = data1['X'].iloc[0]
            ax1.loglog(data1.index, data1[species]/sites,'k--', linewidth=2.5)
        sites = data2['X'].iloc[0]
        ax1.loglog(data2.index, data2[species]/sites,'-')
        plt.xlim(1e-6,1)
        plt.ylim(ymin=1e-6)
        ax1.set_ylabel('site fraction')
        ax1.set_xlabel('time (s)')
        ## Manually constructed label
        x = label_positions[species] if label_positions and species in label_positions else 3.0
        x = 10**(-1*x)
        y = data2[x:][species].iloc[0] / sites
        if x==1: x=1.3 # move them right just a little
        ax1.text(x,y,label, fontsize=10,
                 verticalalignment='center' if x==1.3 else 'bottom',
                 color=sns.color_palette()[i])  

label_positions = {'HX':0.5, 'COX':0, 'CX':2, 'C2HOX':0,
                  'CH3CXO':0, 'CH3CX':0, 'OX':5, 'CHO2X':0,
                  'CH2X':2.5}
fig = plt.figure(figsize=(5,4))
with sns.color_palette("Dark2", len(surface_phase)):
    ax1 = multi_comparison_loglogplot(species_list=surface_phase, label_positions=label_positions)
plt.tight_layout()
plt.savefig('Surface comparison loglog.pdf',bbox_inches='tight')


# # Effect of tolerance
# Now we investigate the effect of gradually decreasing (tightening) the tolerance

# In[30]:

# First, make a series of input files in separate directories

base_directory = 'ch4_co2_tolerances'

with open(os.path.join(base_directory, 'input.template.py')) as infile:
    input_file = infile.read()

def directory(i):
    return os.path.join(base_directory, "ch4_co2_tolerance_m{}".format(i))

for i in range(9):
    tolerance = 0.1**i
    tolerance_string = 'toleranceMoveToCore={:.1e},'.format(tolerance)
    print tolerance_string
    input_file = re.sub('toleranceMoveToCore\s*=\s*(.*?),', tolerance_string, input_file)
    os.path.exists(directory(i)) or  os.makedirs(directory(i))
    with open(os.path.join(directory(i), 'input.py'), 'w') as outfile:
        outfile.write(input_file)
    print "Saved to {}/input.py".format(directory(i))


# In[31]:

# Now run all the jobs
# Don't execute this cell unless you have a while to wait.
import subprocess
import sys
for i in range(9):
    print "Attempting to run job {} in directory {}".format(i, directory(i))
    try:
        retcode = subprocess.call("python $RMGpy/rmg.py {}/input.py".format(directory(i)), shell=True)
        if retcode < 0:
            print >>sys.stderr, "Process was terminated by signal", -retcode
        elif retcode > 0:
            print >>sys.stderr, "Process returned", retcode
        else:
            print "Success"
    except OSError as e:
        print >>sys.stderr, "Execution failed:", e


# In[32]:

# Now read the ends of the log files and extract the mechanism sizes.
epsilon = []
core_species = []
core_rxns = []
edge_species = []
edge_rxns = []

for i in range(9):
    dirname = directory(i)
    print "reading from ", directory(i)
    last_lines = subprocess.check_output(['tail', '-n','6', os.path.join(directory(i),'RMG.log')])

    match = re.search('The final model core has (\d+) species and (\d+) reactions', last_lines)
    if match is None: 
        print "Trouble with {}/RMG.log:\n{}".format(directory(i),last_lines)
    core_species.append(int(match.group(1)))
    core_rxns.append(int(match.group(2)))
    match = re.search('The final model edge has (\d+) species and (\d+) reactions', last_lines)
    if match is None: 
        print "Trouble with {}/RMG.log:\n{}".format(directory(i),last_lines)
    edge_species.append(int(match.group(1)))
    edge_rxns.append(int(match.group(2)))
    epsilon.append( 0.1**i )

#remove the four inerts, N2, Ar, Ne, and He that RMG automatically adds
core_species = np.array(core_species) - 4
edge_species = np.array(edge_species) - 4


# In[33]:

# Now count the reactions by type
Deutschmann = []
abstraction = []
adsorption = []
dissociation = []
for i in range(9):
    ckfilepath = os.path.join(directory(i), 'chemkin', 'chem_annotated.inp')
    re_library = re.compile('! Library reaction: (.*)')
    re_family = re.compile('! Template reaction: (.*)')
    from collections import Counter
    counts = Counter()
    with open(ckfilepath) as ckfile:
        for line in ckfile:
            match = re_family.match(line) or re_library.match(line)
            if match is None:
                continue
            source = match.group(1)
            counts[source] += 1
    counts

    Deutschmann.append(counts['Deutschmann_Ni'])
    abstraction.append(counts['Surface_Abstraction'])
    adsorption.append(counts['Surface_Adsorption_Dissociative'] + counts['Surface_Adsorption_Single'])
    dissociation.append(counts['Surface_Dissociation'])
print 'Deutschmann',Deutschmann
print 'abstraction',abstraction
print 'adsorption',adsorption
print 'dissociation',dissociation

Deutschmann = np.array(Deutschmann)
abstraction = np.array(abstraction)
adsorption = np.array(adsorption)
dissociation = np.array(dissociation)
total = Deutschmann + abstraction + adsorption + dissociation
assert (total == np.array(core_rxns)).all(), "Sum of counters doesn't equal core_rxns from above"


# In[34]:

# Now plot the figure
fig = plt.figure()
fig.set_size_inches(7,3.5) #set size, in inches
gs = matplotlib.gridspec.GridSpec(1, 3)
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1])
ax2 = plt.subplot(gs[2])

ax0.loglog(epsilon, edge_species, 'r^-', label='edge')
ax0.loglog(epsilon, core_species, 'bo-', label='core')

ax1.loglog(epsilon, edge_rxns, 'r^-', label='edge')
ax1.loglog(epsilon, core_rxns, 'bo-', label='core')

#ax2.semilogx(epsilon, Deutschmann, 'k:', label='Deutschmann')
#ax2.semilogx(epsilon, abstraction, 'b-', label='abstraction')
#ax2.semilogx(epsilon, dissociation, 'g--', label='dissociation')
#ax2.semilogx(epsilon, adsorption, 'r-.', label='adsorption')

stacks = ax2.stackplot(epsilon, Deutschmann, adsorption, dissociation, abstraction,
             labels='Deutcschmann Adsorption Dissociation  Abstraction'.split(),
             colors=sns.color_palette("Set2",4)
             )
hatches = ['', '///', '|||', '---']
for i, patch in enumerate(stacks):
    patch.set_hatch(hatches[i])
ax2.set_xscale('log')


#ax0.set_ylim([20,1000])
#ax1.set_ylim([30,2000])
#ax2.set_ylim([0,300])

ax0.set_xlabel('error tolerance, $\epsilon$', fontsize=11)
ax1.set_xlabel('error tolerance, $\epsilon$', fontsize=11)
ax2.set_xlabel('error tolerance, $\epsilon$', fontsize=11)

#ax0.set_ylabel('species')
#ax1.set_ylabel('reactions')
#ax2.set_ylabel('reactions in core')

ax0.set_title('no. of species', fontsize=11)
ax1.set_title('no. of reactions', fontsize=11)
ax2.set_title('no. of core reactions', fontsize=11)

ax0.legend(loc='upper left', bbox_to_anchor=[0.0, 0.9], fontsize=10, frameon=False)
ax1.legend(loc='upper left', bbox_to_anchor=[0.0, 0.9], fontsize=10, frameon=False)
# For the stacked plot we want to plot the legend in reverse order 
# so the bottom patch is on the bottom of the legend
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles[::-1], labels[::-1],
          loc='upper left', bbox_to_anchor=[0.0, 0.9], fontsize=10, frameon=False)


ax0.set_xticks([1.0E-8, 1.0E-4,  1.0])
ax1.set_xticks([1.0E-8, 1.0E-4,  1.0])
ax2.set_xticks([1.0E-8, 1.0E-4,  1.0])
ax0.tick_params(labelsize=10)
ax1.tick_params(labelsize=10)
ax2.tick_params(labelsize=10)

ax0.invert_xaxis()
ax1.invert_xaxis()
ax2.invert_xaxis()

ax0.text(0.1, 0.9, '(a)', fontsize=10,
verticalalignment='bottom', horizontalalignment='left',
transform=ax0.transAxes)
ax1.text(0.1, 0.9, '(b)', fontsize=10,
verticalalignment='bottom', horizontalalignment='left',
transform=ax1.transAxes)
ax2.text(0.1, 0.9, '(c)', fontsize=10,
verticalalignment='bottom', horizontalalignment='left',
transform=ax2.transAxes)

fig.tight_layout()
fig.savefig('mechanism_size.pdf', bbox_inches='tight')


# In[ ]:



