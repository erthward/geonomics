#!/usr/bin/python
#stats.py


'''
Classes and functons to implement calculation and output of statistics
'''

#geonomics imports
from geonomics.utils.io import (_append_array2d_to_array_stack,
                                _append_row_to_csv, _write_dict_to_csv)
from geonomics.ops.selection import _calc_fitness
from geonomics.utils.viz import _check_display

#other imports
import numpy as np
from scipy.stats.stats import pearsonr
from collections import Counter as C
import os
import matplotlib as mpl
_check_display()
import matplotlib.pyplot as plt


######################################
# -----------------------------------#
# CLASSES ---------------------------#
# -----------------------------------#
######################################

#a StatsCollector class, to parameterize and manage calculation
#and collection of stats, then write them to file at the end of 
#each model iteration
class _StatsCollector:
    def __init__(self, model_name, params):

        #set model_name
        self.model_name = model_name

        #set total model time
        self.T = params.model.T

        #grab the stats parameters
        stats_params = params.model.stats

        #a dictionary to link the stats' names in the params dict 
        #to the functions to be called to calculate them 
        self.calc_fn_dict = {'Nt': _calc_Nt,
                             'ld':  _calc_ld,
                             'het': _calc_het,
                             'maf': _calc_maf,
                             'mean_fit': _calc_mean_fitness,
                              }

        #a dictionary to link the stats' names in the params dict 
        #to the functions to be called to write them to disk
        self.write_fn_dict = {'ld':  self._write_array_to_stack,
                             'het': self._write_row_to_csv,
                             'maf': self._write_row_to_csv,
                             }

        #a dictionary to link stats to the file extensions that 
        #should be used to write them to disk
        self.file_suffix_dict = {'Nt': 'OTHER_STATS.csv',
                                'ld':  'LD.txt',
                                'het': 'HET.csv',
                                'maf': 'MAF.csv',
                                'mean_fit': 'OTHER_STATS.csv',
                                 }

        #get the species names
        spps_with_wout_genomes = {str(k):('gen_arch' in v.keys()) for k, v
                                 in params.comm.species.items()}

        #list stats that cannot be calculated for species without genomes
        stats_invalid_wout_genomes = ['ld', 'het', 'maf', 'mean_fit']

        #create a stats attribute, to store all stats calculated
        self.stats = {}
        for spp_name, genome in spps_with_wout_genomes.items():
            self.stats[spp_name] = {}
            for stat, stat_params in stats_params.items():
                #skip species without genomes for stats that need genomes
                if not genome and stat in stats_invalid_wout_genomes:
                    break
                #each spp gets a subdict
                else:
                    #each subdict gets a key for each stat to be calculated
                    if stat_params.calc:
                        #create a subdictionary for each stat, with a list of
                        #NaNs self.T items long, which will be filled in for
                        #each whenever it is sampled (NOTE: this forces all
                        #stats to have the same length so that they all fit
                        #into one pd.DataFrame at the end, and so that plots
                        #easily line up on the same timeframe)
                        self.stats[spp_name][stat]= {
                            'vals': [np.nan]*self.T,
                            'freq': stat_params.freq,
                            #add a 'filepath' key, whose value will be updated
                            #to contain to correct filepaths for each stat
                            'filepath': None,
                            #create tuple of other, stat-specific parameters, 
                            #to later be unpacked as arguments to
                            #the appropriate stat function
                            'other_params': dict([(k,v) for k,v in
                              stat_params.items() if k not in ['calc',
                                                                    'freq']])
                            }

                        #if the freq value is 0, change it to self.T -1, so
                        #that it collects only on the first and last timestep
                        if self.stats[spp_name][stat]['freq'] == 0:
                            self.stats[spp_name][stat]['freq'] = self.T-1

    #create a master method, to be called each timestep, which will make a list 
    #of all stats that need to be calculated that timestep (based on the 
    #calculation-frequencies provided in the params dicts), and then calls the
    #functions to calculate them all and adds the results to self.stats
    def _calc_stats(self, community, t, iteration):
        #set the filepaths, if this is the first timestep of the model
        #iteration
        if t == 0:
            self._set_filepaths(iteration)
        #for each species
        for spp in community.values():
            #list the stats to be calculated this timestep
            if t == self.T-1:
                #calculate all, if it's the last timestep
                calc_list = [*self.stats[spp.name]]
            else:
                #or else only calculate based on the parameterized frequencies
                #for each stat
                calc_list = [k for k,v in self.stats[spp.name].items() if (
                                                        t % v['freq'] == 0)]
            #then calculate each stat
            for stat in calc_list:
                vals = self.calc_fn_dict[stat](spp,
                                **self.stats[spp.name][stat]['other_params'])
                #and add each stat to the right location (by timestep)
                #in its list
                try:
                    self.stats[spp.name][stat]['vals'][t] = vals
                #unless the list isn't long enough (which happens if mod.walk
                #has been used to run the model past its initially stipulated
                #length of time), in which case make it long enough and make
                #the last value the stat just calculated
                except IndexError:
                    stats_list = self.stats[spp.name][stat]['vals']
                    stats_list.extend([np.nan] * (t-len(stats_list)) + [vals])
        #and write whichever stats are necessary to file
        self._write_stats(t)

    #a method to make the filenames for all of the stats to be saved
    def _set_filepaths(self, iteration):
        #get the directory name for this model and iteration
        dirname = os.path.join('GNX_mod-%s' % self.model_name,
                               'it-%i' % iteration)
        #for each species
        for spp_name in [*self.stats]:
            #get the subdirectory name and filename for this species
            subdirname = os.path.join(dirname, 'spp-%s' % spp_name)
            #make this subdir, and any parent dirs as necessary
            os.makedirs(subdirname, exist_ok = True)
            #create the filename and filepath for this spp, for each stat
            for stat in [*self.stats[spp_name]]:
                filename = 'mod-%s_it-%i_spp-%s_%s' % (self.model_name,
                    iteration, spp_name, self.file_suffix_dict[stat])
                filepath = os.path.join(subdirname, filename)
                #add the filepath for this stat to self.stats
                self.stats[spp_name][stat]['filepath'] = filepath

    #wrapper around io.append_array2d_to_array_stack
    #TODO WHAT TO DO WITH t IN THIS CASE?? CAN'T ADD TO txt 3D ARRAY FILE
    def _write_array_to_stack(self, filepath, array, t):
        _append_array2d_to_array_stack(filepath, array)

    #wrapper around io.append_row_to_csv
    def _write_row_to_csv(self, filepath, array, t):
        _append_row_to_csv(filepath, array, t)

    #use io._write_dict_to_csv to write to disk all "other stats", i.e. 
    #all stats that collect only a single value per species per timestep
    #TODO: CHANGE THE 'OTHER STATS' NAMING CONVENTION TO SOMETING MORE
    #DESCRIPTIVE
    def _write_other_stats(self):
        for spp, spp_stats in self.stats.items():
            #get a dictionary of the data values for all stats that are to be
            #written just once at the end of the iteration
            data_dict = {k:v['vals'] for k,v in spp_stats.items() if
                              'OTHER_STATS' in v['filepath']}
            #they all have the same filepath, so just grab the first
            filepath = [*spp_stats.values()][0]['filepath']
            #write to disk
            _write_dict_to_csv(filepath, data_dict)

    #method to write stats to files, in the appropriate directory (by model
    #and iteration number), and with the appropriate spp names in the filenames
    def _write_stats(self, t):
        #for each species
        for spp_name, spp_stats in self.stats.items():
            #for each stat
            write_list = [k for k,v in spp_stats.items() if t % v['freq'] == 0]
            for stat, stat_dict in spp_stats.items():
                #get the filepath
                filepath = stat_dict['filepath']
                #if the filepath does not contain "OTHER_STATS" then it is a
                #stat that produces more than a single value per species per
                #timestep it is collected, so write the data to disk
                #intermittently and then delete the data from memory (if it was
                #collected this timestep)
                if "OTHER_STATS" not in filepath and stat in write_list:
                    #get the correct write_fn for this stat
                    write_fn = self.write_fn_dict[stat]
                    #call the write_fn to write the data to disk
                    write_fn(filepath, stat_dict['vals'][t], t)
                    #then replace the last data collected prior to this
                    #timestep's data with None, to free up memory but still
                    #maintain the latest data in case of plotting
                    rev_nonnull = [n for n, v in enumerate(
                        stat_dict['vals'][::-1]) if (v is not np.nan and
                                                        v is not None)]
                    nonnull = [range(len(
                        stat_dict['vals']))[::-1][n] for n in rev_nonnull]
                    nonnull = [v for v in nonnull if v != t]
                    for v in nonnull:
                        stat_dict['vals'][v] = None


        #or write all 'other stats' to disk, if it's the last timestep
        if t == self.T-1:
            self._write_other_stats()

    #method to plot whichever stat as a function of runtime
    def _plot_stat(self, stat, spp_name=None):
        #check that the stat argument is valid
        assert type(stat) is str, "The 'stat' argument must be a string."
        assert stat in [*self.stats.values()][0].keys(), ("The 'stat' "
            "argument must name a statistic that was collected. Valid values: "
            "%s.") % (','.join(["'%s'" % val for val in
                                [*self.stats.values()][0].keys()]))
        #get the list of spps to plot
        if spp_name is None:
            spp_names = [*self.stats]
        elif (spp_name is not None
              and type(spp_name) is str and spp_name in [*self.stats]):
            spp_names = [spp_name]
        else:
            raise ValueError(("The 'spp_name' argument, if provided, "
                "must be a string containing a valid species name."))
        #create the figure
        fig = plt.figure()
        #plot each species for the chosen statistic 
        for n, spp_name in enumerate(spp_names):
            #get the stat values to plot
            vals = self.stats[spp_name][stat]['vals']
            #plot 'Nt' or 'mean_fit'
            if stat in ['Nt', 'mean_fit']:
                #add axes objects horizontally across
                ax = fig.add_subplot(1, len(spp_names), n+1)
                #get the indices of non-NaN values to be plotted
                indices_to_plot = np.array(np.where(
                                    np.invert(np.isnan(vals)))[0])
                #get the timesteps at the non-NaN values
                x = np.arange(0, len(vals))[indices_to_plot]
                #get the non-NaN values
                y = np.array(vals)[indices_to_plot]
                #plot a dotted line (which necessarily linearly interpolates 
                #between collected timesteps if not all timesteps
                #were collected)
                plt.plot(x, y, ':')
                #and plot dots at each of the collected timesteps
                plt.plot(x, y, '.')
                #set the title to the stat and the species' name
                ax.set_title("SPP: '%s'" % (spp_name))
                #set the x- and y-labels
                plt.xlabel('timestep')
                plt.ylabel(stat)

            #or plot 'maf' or 'het'
            elif stat in ['het', 'maf']:
                #add axes objects horizontally across
                ax = fig.add_subplot(1, len(spp_names), n+1)
                #get the reversed-list index of the last set of values 
                #calculated
                rev_idx_last_vals = [n for n,v in enumerate(vals[::-1]) if (
                                    v is not None and v is not np.nan)][0]
                #get the last set of values calculated
                last_vals = vals[::-1][rev_idx_last_vals]
                #get the timestep of the last set of values
                t_last_vals = range(len(vals))[::-1][rev_idx_last_vals]
                #plot the values
                plt.plot(range(len(last_vals)), last_vals, '-')
                #set the title to the species' name and timestep of the
                #values plotted
                ax.set_title("SPP: '%s';   T: %i" % (spp_name, t_last_vals))
                #set the x- and y-labels
                plt.xlabel('locus')
                plt.ylabel(stat)

            #or plot 'ld'
            elif stat in ['ld']:
                #get the reversed-list index of the last set of values 
                #calculated
                rev_idx_last_vals = [n for n,v in enumerate(vals[::-1]) if (
                    v is not None and v is not np.nan)][0]
                #get the last set of values (i.e. r^2 array) calculated
                r2_mat = vals[::-1][rev_idx_last_vals]
                #get the timestep of the last set of values
                t_last_vals = range(len(vals))[::-1][rev_idx_last_vals]
                #add axes objects horizontally across, in two rows
                ax = fig.add_subplot(2, len(spp_names), n+1)
                #plot the LD matrix in row 1
                plt.imshow(np.clip(r2_mat, a_min = 0, a_max = None),
                                            interpolation = 'nearest')
                plt.colorbar()
                #set plot title
                ax.set_title(("SPP: '%s';   T: %i\nLocus-wise "
                            "linkage matrix") %  (spp_name, t_last_vals))
                #set the x- and y-labels
                plt.xlabel('locus')
                plt.ylabel('locus')
                ax = fig.add_subplot(2, len(spp_names), n+1+len(spp_names))
                #plot of mean linkage values
                r2_list = [r2_mat[0,1]]
                L = r2_mat.shape[0]
                for i in range(1,L-1):
                    r2_list.append(np.mean([r2_mat[i-1,i], r2_mat[i,i+1]]))
                r2_list.append(r2_mat[L-2,L-1])
                plt.scatter(range(L), r2_list, c = 'red', marker = 'o', s=25)
                #set plot title
                ax.set_title("Locus-wise mean linkage values")
                #set the x- and y-labels
                plt.xlabel('locus')
                plt.ylabel('mean linkage')

            #or else return informative error message
            else:
                raise ValueError(("The value provided for the 'stat' argument "
                    "is not a valid statistic. Valid values include: %s\n\n")%(
                    ','.join(['%s' % k for k in [*self.calc_fn_dict]])))
        #set the main title to the stat plotted
        fig.suptitle('STAT: %s' % stat)
        #show the image
        fig.show()


######################################
# -----------------------------------#
# FUNCTIONS -------------------------#
# -----------------------------------#
######################################

#method to get pop size (NOTE: not actually calculating it)
def _calc_Nt(spp):
    Nt = spp.Nt[-1]
    return(Nt)


def _calc_ld(spp, plot = False):

    #TODO: I should also include (either as an alternative within this fn,
    #or as separate fn) the option to calculate D'

    #TODO: I keep getting warnings like the following, which could just be 
    #due to divison of small floating-point numbers, but I should figure out 
    #exactly what's going on and be sure everything checks out. WARNING:
    # stats.py:117: RuntimeWarning: invalid value encountered in double_scalars

    speciome = spp._get_genotypes()
    n = np.shape(speciome)[0] #num individs
    x = np.shape(speciome)[2] #ploidy
    N = n*x
    L = spp.gen_arch.L
    assert L == np.shape(speciome)[1], ("The length of the 1st dimension "
                            "of speciome doesn't equal spp.genomic_arch.L")

    r2_mat = np.zeros([L]*2) * np.nan # vals default to NaN

    for i in range(L):
        for j in range(i+1, L):
            #calculate freq of allele 1 at locus i
            f1_i = np.sum(speciome[:,i,:], axis = None)/(N)
            #calculate freq of allele 1 at locus j
            f1_j = np.sum(speciome[:,j,:], axis = None)/(N)
            #calculate freq of chroms with 1_1 haplotype at loci i and j
            f11_ij = float(np.sum(speciome[:,[i,j],:].sum(axis = 1) ==2,
                                                        axis = None))/(N)
            D_1_1 = f11_ij - (f1_i * f1_j)
            r2 = (D_1_1**2)/(f1_i*(1-f1_i)*f1_j*(1-f1_j))
            # add to both triangular halves, to produce a symmetric matrix
            r2_mat[i,j] = r2
            r2_mat[j,i] = r2

    return(r2_mat)


#function to calculate the locus-wise (if mean == False) or mean (if
#mean == True) heterozygosity of the species
def _calc_het(spp, mean=False):
    #get pop size
    N = len(spp)
    #get the speciome
    speciome = spp._get_genotypes()
    #calculate the frequency of heterozygotes, locus-wise
    het = np.sum(np.mean(speciome, axis = 2) == 0.5, axis = 0)/N
    #get the mean heterozygosity, if mean argument is True
    if mean:
        het = mean(het)
    return(het)

#function to calculate the locus-wise minor allele frequency of the species
def _calc_maf(spp):
    #get two times the pop size
    two_N = 2*len(spp)
    #get the speciome
    speciome = spp._get_genotypes()
    #get the frequencies of 1-alleles for all loci
    freqs_1 = np.sum(np.sum(speciome, axis = 2), axis = 0)/two_N
    #find all loci where the 1-allele is the major allele
    majors = np.where(freqs_1 > 0.5)
    #replace the locations where 1 is the major allele with 0-allele freq
    maf = freqs_1[:]
    maf[majors] = 1 - freqs_1[majors]
    return(maf)


#function to calculate the mean fitness of the species
def _calc_mean_fitness(spp):
    #calculate the mean fitness, if this species has traits
    if spp.gen_arch.traits is not None:
        mean_fit = np.mean(_calc_fitness(spp))
    #or else return NaN
    else:
        mean_fit = np.nan
    return(mean_fit)
