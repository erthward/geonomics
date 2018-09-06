#!/usr/bin/python
#stats.py


'''
##########################################

Module name:              stats

Module contents:          - definition of Stats class (i.e. structured container
                            for stats calculated during model run)
                          - definition of functions for calculating various stats,
                            at specified frequencies and with specified arguments,
                            according to the contents of the params.model.stats section


Author:                   Drew Ellison Hart
Email:                    drew.hart@berkeley.edu
Github:                   URL
Start date:               01-01-18
Documentation:            URL


##########################################
'''

#geonomics imports


#other imports
import numpy as np
from scipy.stats.stats import pearsonr
from collections import Counter as C
import matplotlib.pyplot as plt


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

#a StatsCollector class, to parameterize and manage calculation
#and collection of stats, then write them to file at the end of 
#each model iteration
class StatsCollector:
    def __init__(self, model_name, params):

        #set model_name
        self.model_name = model_name

        #set total model time
        self.T = params.model.time.T

        #grab the stats parameters
        stats_params = params.model.stats

        #create a dictionary to link the stats' names in the params dict 
        #to the functions to be called by them (defined in stats.py)
        self.calc_fn_dict = {  'Nt': calc_Nt,
                                'ld':  calc_ld,
                                'het': calc_het,
                                'maf': calc_maf,
                                'mean_fit': calc_mean_fitness,
                              }

        #get the population names
        pop_names = [v.init.name for v in params.comm.pops.values()]

        #list stats that cannot be calculated for populations without genomes
        stats_invalid_wout_genomes = ['ld', 'het', 'maf', 'mean_fit']

        #create a stats attribute, to store all stats calculated
        self.stats = {}
        for pop_name in pop_names:
            #skip populations without genomes for stats that require genomes
            if pop.gen_arch is None and stat in stats_invalid_wout_genomes:
                break
            #each pop gets a subdict
            else:
                self.stats[pop_name] = {}
                for stat, stat_params in stats_params.items():
                    #each subdict gets a key for each stat to be calculated
                    if stat_params.calc:
                        #create a subdictionary for each stat, with a list of NaNs
                        #self.T items long, which will be filled in for each 
                        #whenever it is sampled (NOTE: this forces all stats to 
                        #have the same length so that they all fit into one 
                        #pd.DataFrame at the end, and so that plots easily line 
                        #up on the same timeframe)
                        self.stats[pop_name][stat]= {
                            'values': [np.nan]*self.T,
                            #'values': np.float32([np.nan]*self.T),
                            'freq': stat_params.freq,
                            #create tuple of other, stat-specific parameters, 
                            #to later be unpacked as arguments to the appropriate 
                            #stat function
                            'other_params': dict([(k,v) for k,v in
                                    stat_params.items() if k not in ['calc', 'freq']])
                            }

                        #if the freq value is 0, change it to self.T -1, so that it
                        #collects only on the first and last timestep
                        if self.stats[pop_name][stat]['freq'] == 0:
                            self.stats[pop_name][stat]['freq'] = self.T-1

    #create a master method, to be called each timestep, which will make a list 
    #of all stats that need to be calculated that timestep (based on the 
    #calculation-frequencies provided in the params dicts), and then calls the
    #functions to calculate them all and adds the results to self.stats
    def calc_stats(self, community, t, iteration):
        for pop in community.values():
            #list the stats to be calculated this timestep
            calc_list = [k for k,v in self.stats[pop.name].items() if t%v['freq'] == 0]
            #then calculate each and append it to self.stats
            for stat in calc_list:
                vals = self.calc_fn_dict[stat](pop,
                                        **self.stats[pop.name][stat]['other_params'])
                self.stats[pop.name][stat]['values'][t] = vals
        #and write the stats to file, if t = self.T-1
        if t == self.T-1:
            self.write_stats(iteration)


    #TODO: THIS OF COURSE WON'T WORK BECAUSE MOST STATS ARE NOT A SINGLE VALUE
    #PER POP PER TIMESTEP; NEED TO RETHINK
    #method to write stats to files, in the appropriate directory (by model
    #and iteration number), and with the appropriate pop names in the filenames
    def write_stats(self, iteration):
        #get the directory name for this model and iteration
        dirname = os.path.join('GEONOMICS_mod-%s' % self.model_name, 
                               'it-%i' % iteration)
        for pop_name in [*self.stats]:
            #get the subdirectory name and filename for this population
            subdirname = os.path.join(dirname, 'pop-%s' % pop_name)
            filename = 'mod-%s_it-%i_pop-%s_STATS.csv'
            filepath = os.path.join(subdirname, filename)
            #make a pandas DataFrame from the population's stats
            df = pd.DataFrame.from_dict({k:v['values'] for k,v in
                                         self.stats[pop_name].items()})
            #add a timestep column
            df['t'] = [*df.index]
            #reorder the columns so that timestep is the first
            ordered_cols = ['t'] + [col for col in df.columns if col != 't']
            df = df[ordered_cols]
            #write to disk
            df.to_csv(path_or_buf= filepath, float_format= '0.5f', index= False)


    #TODO: USE A CHECK OF THE SHAPES OF THE STATS TO DECIDE HOW TO PLOT
    #method to plot whichever stat as a function of runtime
    def show_stat(self, stat, pop_name=None):
        #check that the stat argument is valid
        assert type(stat) is str, "The 'stat' argument must be a string."
        assert stat in [*self.stats.values()][0].keys(), ("The 'stat' "
            "argument must name a statistic that was collected. Valid values: "
            "%s.") % (','.join(["'%s'" % val for val in
                                [*self.stats.values()][0].keys()]))
        #get the list of pops to plot
        if pop_name is None:
            pop_names = [*self.stats]
        elif (pop_name is not None 
              and type(pop_name) is str and pop_name in [*self.stats]):
            pop_names = [pop_name]
        else:
            raise ValueError(("The 'pop_name' argument, if provided, "
                "must be a string containing a valid population name."))
        #create the figure
        fig = plt.figure()
        #plot each population for the chosen statistic 
        for n, pop_name in enumerate(pop_names):
            #add axes objects horizontally across
            ax = fig.add_subplot(1, len(pop_names), n)
            #get the stat values to plot
            vals = self.stats[pop_name]['values']
            #get the indices of non-NaN values to be plotted
            indices_to_plot = np.where(np.invert(np.isnan(vals)))
            #get the timesteps at the non-NaN values
            x = np.arange(0, len(vals))[indices_to_plot]
            #get the non-NaN values
            y = vals[indices_to_plot]
            #plot a dotted line (which necessarily linearly interpolates 
            #between collected timesteps if not all timesteps were collected)
            plt.plot(x, y, ':')
            #and plot dots at each of the collected timesteps
            plt.plot(x, y, '.')
            #set the title to the population's name
            ax.set_title("POP: '%s'" % pop_name)
        #show
        fig.show()


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

#TODO: either get rid of this, or just make it check the pop.Nt list instead
def calc_Nt(pop):
    Nt = len(pop)
    return(Nt)


#TODO: REVIEW AND FIX/EDIT AS NEEDED
#TODO: INCORPORATE THE PLOTTING STUFF BELOW INTO self.show_stat() INSTEAD
def calc_ld(pop, plot = False):
    
    #TODO: I should also include (either as an alternative within this fn,
    #or as separate fn) the option to calculate D'

    #TODO: I keep getting warnings like the following, which could just be 
    #due to divison of small floating-point numbers, but I should figure out 
    #exactly what's going on and be sure everything checks out. WARNING:
    # stats.py:117: RuntimeWarning: invalid value encountered in double_scalars

    populome = get_populome(pop)
    #populome = np.array([ind.genome for ind in pop.values()])
    n = np.shape(populome)[0] #num individs
    x = np.shape(populome)[2] #ploidy
    N = n*x
    L = pop.gen_arch.L
    assert L == np.shape(populome)[1], ("The length of the 1th dimension "
                            "of populome doesn't equal pop.genomic_arch.L")

    r2_mat = np.zeros([L]*2)-1 # -1 serves as a 'no data' value here

    for i in range(L):
        for j in range(i+1, L):
            #calculate freq of allele 1 at locus i
            f1_i = np.sum(populome[:,i,:], axis = None)/(N)
            #calculate freq of allele 1 at locus j
            f1_j = np.sum(populome[:,j,:], axis = None)/(N)
            #calculate freq of chroms with 1_1 haplotype at loci i and j
            f11_ij = float(np.sum(populome[:,[i,j],:].sum(axis = 1) ==2,
                                                        axis = None))/(N)
            D_1_1 = f11_ij - (f1_i * f1_j)
            r2 = (D_1_1**2)/(f1_i*(1-f1_i)*f1_j*(1-f1_j))
            r2_mat[i,j] = r2

    if plot == True:
        fig = plt.figure()
        ax = fig.add_subplot(1,2,1)
        #plot of LD matrix
        plt.imshow(np.clip(r2_mat, a_min = 0, a_max = None),
                                    interpolation = 'nearest')
        ax = fig.add_subplot(1,2,2)
        #plot of mean linkage values
        r2_list = [r2_mat[0,1]]
        for i in range(1,L-1):
            r2_list.append(np.mean([r2_mat[i-1,i], r2_mat[i,i+1]]))
        r2_list.append(r2_mat[L-2,L-1])
        plt.scatter(range(L), r2_list, c = 'red', marker = 'o', s=25)
    return(r2_mat)


#function to calculate the locus-wise (if mean == False) or mean (if
#mean == True) heterozygosity of the population
def calc_het(pop, mean=False):
    #get pop size
    N = len(pop)
    #get the populome
    populome = get_populome(pop)
    #calculate the frequency of heterozygotes, locus-wise
    het = np.sum(np.mean(populome, axis = 2) == 0.5, axis = 0)/N
    #get the mean heterozygosity, if mean argument is True
    if mean:
        het = mean(het)
    return(het)

#function to calculate the locus-wise minor allele frequency of the population
def calc_maf(pop):
    #get two times the pop size
    two_N = 2*len(pop)
    #get the populome
    populome = get_populome(pop)
    #get the frequencies of 1-alleles for all loci
    freqs_1 = np.sum(np.sum(populome, axis = 2), axis = 0)/two_N
    #find all loci where the 1-allele is the major allele
    majors = np.where(freqs_1 > 0.5)
    #replace the locations where 1 is the major allele with 0-allele freq
    maf = freqs_1[:]
    maf[majors] = 1 - freqs_1[majors]
    return(maf)

#function to calculate the mean fitness of the population
def calc_mean_fitness(pop):
    pass

#helper function for creating a 3-d array of all individuals' genotypes
#(a 'populome')
def get_populome(pop):
    populome = np.stack([ind.genome for ind in pop.values()])
    return(populome)


