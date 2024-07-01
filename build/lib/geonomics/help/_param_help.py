#!/usr/bin/python
# _param_help.py

# flake8: noqa


'''
Helper functions for assessing the effects of different parameter settings.
'''

# geonomics imports
from geonomics.utils.viz import _check_display

# other imports
import os
import numpy as np
import matplotlib as mpl
_check_display()
import matplotlib.pyplot as plt
from scipy.stats import levy as _s_levy
from copy import deepcopy

# ------------------------------------
# CLASSES ---------------------------
# ------------------------------------


# --------------------------------------
# FUNCTIONS ---------------------------
# --------------------------------------


# functions to visualize the different distributions that can be
# parameterized from the parameters file
def plot_distr_movement_distance(spp=None, distance_distr_param1=None,
                                 distance_distr_param2=None,
                                 distance_distr='levy'):
    if spp is not None:
        distance_distr_param1 = spp.movement_distance_distr_param1
        distance_distr_param2 = spp.movement_distance_distr_param2
    else:
        assert (distance_distr_param1 is not None and
                distance_distr_param2 is not None), ('If a Species object '
                                                     'is not provided then '
                                                     'the parameter values '
                                                     'must be provided.')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if distance_distr == 'levy':
        fig.suptitle(('movement distance: '
                  '~Levy($\loc$=%.4E, $\scale$=%.4E)') % (distance_distr_param1,
                                                        distance_distr_param2))
        vals = _s_levy.rvs(loc=distance_distr_param1,
                              scale=distance_distr_param2, size=10000)
    elif distance_distr == 'wald':
        fig.suptitle(('movement distance: '
                  '~Wald($\mu$=%.4E, $\scale$=%.4E)') % (distance_distr_param1,
                                                         distance_distr_param2))
        vals = np.random.wald(mean=distance_distr_param1,
                              scale=distance_distr_param2, size=10000)
    ax.hist(vals, bins = 25)
    ax.set_xlim((min(vals), max(vals)))
    plt.show()


def plot_distr_movement_direction(spp=None, direction_distr_mu=None,
                            direction_distr_kappa=None):
    if spp is not None:
        direction_distr_mu = spp.direction_distr_mu
        direction_distr_sigma = spp.direction_distr_kappa
    else:
        assert (direction_distr_mu is not None and
                direction_distr_kappa is not None), ('If a Species object '
                                                     'is not provided then '
                                                     'the parameter values '
                                                     'must be provided.')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.suptitle(('movement direction: '
                 '~von Mises($\mu$ = %.4E, '
                  '$\kappa$ = %.4E)') % (direction_distr_mu,
                                         direction_distr_kappa))
    vals = np.random.vonmises(direction_distr_mu, direction_distr_kappa, 10000)
    ax.hist(vals, bins = 25)
    ax.set_xlim((min(vals), max(vals)))
    plt.show()


def plot_distr_dispersal_distance(spp=None, dispersal_distr_param1=None,
                                  dispersal_distr_param2=None,
                                  dispersal_distance_distr='levy'):
    if spp is not None:
        dispersal_distr_param1 = spp.dispersal_distance_distr_param1
        dispersal_distr_param2 = spp.dispersal_distance_distr_param2
    else:
        assert (dispersal_distr_param1 is not None and
                dispersal_distr_param2 is not None), ('If a Species object '
                                                     'is not provided then '
                                                     'the parameter values '
                                                     'must be provided.')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if dispersal_distance_distr == 'levy':
        fig.suptitle(('dispersal distance: '
                      '~Levy(loc=%.4E, scale=%.4E)') % (dispersal_distr_param1,
                                                           dispersal_distr_param2))
        vals = _s_levy.rvs(dispersal_distr_param1, dispersal_distr_param2, 10000)
    elif dispersal_distance_distr == 'wald':
        fig.suptitle(('dispersal distance: '
                      '~Wald(mean=%.4E, scale=%.4E)') % (dispersal_distr_param1,
                                                           dispersal_distr_param2))
        vals = np.random.wald(dispersal_distr_param1, dispersal_distr_param2, 10000)
    ax.hist(vals, bins = 25)
    ax.set_xlim((min(vals), max(vals)))
    plt.show()


def plot_n_births(spp=None, n_births_distr_lambda=None):
    if spp is not None:
        n_births_distr_lambda = spp.n_births_distr_lambda
    else:
        assert n_births_distr_lambda is not None, ('If a Species object '
                                                  'is not provided then '
                                                  'the parameter values '
                                                  'must be provided.')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.suptitle(('number of births per mating event: '
                  '~Poisson($\lamba$=%.4E)') % n_births_distr_lambda)
    vals = np.random.poisson(n_births_distr_lambda, 10000)
    ax.hist(vals, bins = 25)
    ax.set_xlim((min(vals), max(vals)))
    plt.show()


def plot_distr_delet_effect_sizes(spp=None, delet_alpha_distr_shape=None,
                                  delet_alpha_distr_scale=None):
    if spp is not None:
        delet_alpha_distr_shape = spp.gen_arch.delet_alpha_distr_shape
        delet_alpha_distr_scale = spp.gen_arch.delet_alpha_distr_scale
    else:
        assert (delet_alpha_distr_shape is not None
                and delet_alpha_distr_scale is not None), ('If a Species '
                                                           'object is not '
                                                           'provided then '
                                                           'the parameter '
                                                           'values must be '
                                                           'provided.')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.suptitle(('effect sizes of deleterious mutations: '
                  '~Gamma($shape$=%.4E, '
                  '$scale$=%.4E)') % (delet_alpha_distr_shape,
                                      delet_alpha_distr_scale))
    vals = np.random.gamma(delet_alpha_distr_shape,
                           delet_alpha_distr_scale, 10000)
    ax.hist(vals, bins = 25)
    ax.set_xlim((min(vals), max(vals)))
    plt.show()


def plot_distr_trait_effect_sizes(spp=None, r_distr_alpha=None,
                            r_distr_beta=None):
    if spp is not None:
        r_distr_alpha = spp.gen_arch.recombinations._r_distr_alpha
        r_distr_beta = spp.gen_arch.recombinations._r_distr_beta
    else:
        assert (r_distr_alpha is not None
                and r_distr_beta is not None), ('If a Species '
                                                'object is not '
                                                'provided then '
                                                'the parameter '
                                                'values must be '
                                                'provided.')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.suptitle(('recombination rates: '
                  '~Beta($\alpha$=%.4E, $\beta$=%.4E)') % (r_distr_alpha,
                                                           r_distr_beta))
    vals = np.random.beta(r_distr_alpha, r_distr_beta, 10000)
    ax.hist(vals, bins = 25)
    ax.set_xlim((min(vals), max(vals)))
    plt.show()


def plot_distr_recomb_rates(spp=None, trt_num=None, alpha_distr_mu=None,
                            alpha_distr_sigma=None):
    if spp is not None:
        assert trt_num is not None, ('Trait number must be provided if a '
                                     'Species object is provided.')
        alpha_distr_mu = spp.gen_arch.traits[trt_num]._alpha_distr_mu
        alpha_distr_sigma = spp.gen_arch.traits[trt_num]._alpha_distr_sigma
    else:
        assert (alpha_distr_mu is not None
                and alpha_distr_sigma is not None), ('If a Species '
                                                     'object is not '
                                                     'provided then '
                                                     'the parameter '
                                                     'values must be '
                                                     'provided.')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.suptitle(('effect sizes of loci, trait %i: '
                  '~Normal($\mu$=%.4E, $\sigma$=%.4E)') % (alpha_distr_mu,
                                                           alpha_distr_sigma))
    vals = np.random.normal(alpha_distr_mu, alpha_distr_sigma, 10000)
    ax.hist(vals, bins = 25)
    ax.set_xlim((min(vals), max(vals)))
    plt.show()
