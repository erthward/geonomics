#!/usr/bin/python

#NOTE: run from within ./python_model/scripts_for_practice/profiling


import cProfile, pstats, StringIO, os, sys

#run setup scripts
execfile('../../params.py')
execfile('../PRACTICE_create_land_genomic_arch_pop.py')
execfile('../PRACTICE_imports_and_reloads.py')


#set pop.K
pop.set_K(land.scapes[params['n_movement_surf_scape']])


#define the __profile_main__ function
execfile('./PRACTICE_profile_main.py')

#cProfile.run that function and dump the resulting stats into a file
cProfile.run('__profile_main__()', fileanme = 'profile_test.pstats')

#use gprof2dot to produce a plot of it; https://github.com/jrfonseca/gprof2dot
cmd = 'gprof2dot -f pstats profile_test.pstats | dot -Tpng -o profile_test.png'
os.system(cmd)
