#!/usr/bin/python

import numpy as np
import statsmodels.api as sm
from collections import OrderedDict as OD

def MMRR(Y, X, Xnames=None, nperm=999):
    """
    This is a port of Ian Wang's MMRR script, which lives here:
        https://nature.berkeley.edu/wanglab/data/

    MMRR performs Multiple Matrix Regression with Randomization analysis

    Parameters
    ----------

    Y: {2d numpy array}
       A dependent distance matrix

    X: {[2d numpy array, 2d numpy array, ...]}
        A list of independent distance matrices

    Xnames: {[str, str, ...]}, optional
        A list of variable names for the X matrices (defaults to ["X1", "X2",
        ...])

    nperm: int, optional
       The number of permutations to use for the permutation test
    """
    # get number Y rows, then unfold lower triangular
    nrowsY = Y.shape[0]
    y = _unfold_tril(Y)
    # process X names, then unfold all Xs' lower triangulars
    if Xnames is None:
        Xnames = ["X%i" % i for i in range(1, len(X)+1)]
    xs = [_unfold_tril(x) for x in X]
    # get in correct array formats for OLS
    mod_y, mod_x = _prep_mod_data(y, xs)
    # fit the linear regression and get the stats
    mod = sm.OLS(mod_y, mod_x).fit()
    coeffs = mod.params
    r2 = mod.rsquared
    tstat = mod.tvalues
    Fstat = mod.fvalue
    tprob = np.ones((len(tstat)))
    Fprob = 1
    # get the row numbers
    rownums = [*range(nrowsY)]
    # perform permutations
    for i in range(nperm):
        # shuffle the row numbers
        np.random.shuffle(rownums)
        Yperm = Y[rownums, rownums]
        yperm = _unfold_tril(Yperm)
        permmod_y, permmod_x = _prep_mod_data(yperm, xs)
        permmod = sm.OLS(permmod_y, permmod_x).fit()
        tprob += (permmod.tvalues >= np.abs(tstat))
        Fprob += (permmod.fvalue >= Fstat)

    # calculate the empirical p-values
    tp = tprob/(nperm+1)
    Fp = Fprob/(nperm+1)

    # return values
    coeff_names = ["Intercept"] + Xnames
    output = OD()
    output["R^2"] = r2
    output.update({c: cval for c, cval in zip(coeff_names, coeffs)})
    output.update({c+ "(t)": tval for c, tval in zip(coeff_names, tstat)})
    output.update({c+ "(p)": pval for c, pval in zip(coeff_names, tp)})
    output["F-statistc"] = Fstat
    output["F p-value"] = Fp

    return output


def _prep_mod_data(y, xs, intercept=True):
    """
    put unfolded y and x matrix tril data into correct shapes
    """
    y_out = y.reshape(len(y), 1)
    x_out = np.vstack(xs).T
    if intercept:
        x_out = np.hstack([np.ones(y_out.shape), x_out])
    return y_out, x_out


def _unfold_tril(A):
    """
    unfolds the lower triangular elements of a matrix into a vector (i.e. 1d
    array)
    """
    vec = np.tril(A, k=-1)[np.tril_indices(A.shape[0], k=-1)]
    return vec
