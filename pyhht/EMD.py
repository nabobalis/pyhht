from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from scipy import interpolate
from scipy.stats import pearsonr

__all__ = ['emd']

def exterma(data, extenstion='extrema', n=2):
    """
    Takes an 1D array and finds exterma points.
    """
    N = data.shape[0]
    min_env = np.zeros(N)
    max_env = min_env.copy()
    min_env = np.logical_and(
            np.r_[True, data[1:] < data[:-1]],
            np.r_[data[:-1] < data[1:], True])
    max_env = np.logical_and(
                        np.r_[True, data[1:] > data[:-1]],
                        np.r_[data[:-1] > data[1:], True])
    max_env[-1] = min_env[0] = False
    min_env = min_env.nonzero()[0]
    max_env = max_env.nonzero()[0]

    data_min = data[min_env]
    data_max = data[max_env]

    min_arr = np.array([min_env, data_min])
    max_arr = np.array([max_env, data_max])

    if min_env.shape[0] <= 2 or max_env.shape[0] <= 2:
        #If this IMF has become a straight line
        pass
    elif  extenstion == 'extrema':
        left_min = np.zeros([2,n])
        right_min = np.zeros([2, n])

        left_max = np.zeros([2,n])
        right_max = np.zeros([2, n])

        for i in range(1, n+1):
            left_max[:, i-1] = [-1*min_env[n-i], data_max[n-i]]
            left_min[:, i-1] = [-1*max_env[n-i], data_min[n-i]]

            right_max[:, i-1] = [2*N - min_env[-i], data_max[-i]]
            right_min[:, i-1] = [2*N - max_env[-i], data_min[-i]]

        min_arr = np.concatenate([left_min, min_arr, right_min], axis=1)
        max_arr = np.concatenate([left_max, max_arr, right_max], axis=1)
    else:
        min_arr = np.array([min_env, data_min])
        max_arr = np.array([max_env, data_max])

    return min_arr, max_arr

def envelope(min_arr, max_arr, N, n, periodic=0):
    #Cubic Spline by default
    order_max = 3
    order_min = 3

    min_arr = np.asarray(min_arr)
    max_arr = np.asarray(max_arr)

    if min_arr.shape[1]-n < 4:
        order_min = 1 #Do linear interpolation if not enough points
    elif min_arr.shape[1]-n < 5:
        order_min = 2 #Do quad interpolation if not enough points
    else:
        order_min = 3

    if max_arr.shape[1]-n < 4:
        order_max = 1  #Do linear interpolation if not enough points
    elif max_arr.shape[1]-n < 5:
        order_max = 2 #Do quad interpolation if not enough points
    else:
        order_max = 3
    # Mirror Method requires per flag = 1
    # No extrapolation requires per flag = 0
    t = interpolate.splrep(*min_arr, k=order_min, per=periodic)
    top = interpolate.splev(np.arange(N), t)

    b = interpolate.splrep(*max_arr, k=order_max, per=periodic)
    bot = interpolate.splev(np.arange(N), b)
    mean = (top + bot)/2

    return mean

def emd(data, nimfs=12, extrapolation='exterma', n=2,
        shifting_distance=0.2, pearson=True):
    """
    Perform a Empirical Mode Decomposition on a data set.

    This function will return an array of all the Imperical Mode Functions as
    defined in [1]_, which can be used for further Hilbert Spectral Analysis.

    The EMD uses a spline interpolation function to approcimate the upper and
    lower envelopes of the signal, this routine implements a extrapolation
    routine as described in [2]_ as well as the standard spline routine.
    The extrapolation method removes the artifacts introduced by the spline fit
    at the ends of the data set, by making the dataset a continuious circle.

    Many thousands of papers have been published with ideas to improve the EMD
    procress. One is paper [3]_, that is used for the exterma mirror.

    Parameters
    ----------
    data : array_like
            Signal Data 1-D array.
    extrapolation : str, optional
            Sets the extrapolation method for edge effects.
            Options: None
                     'mirror'
                     'extrema'
            Default: 'mirror'
    n: int, optional
            Sets the number of points used for the exterma mirror method.
    nimfs : int, optional
            Sets the maximum number of IMFs to be found
            Default : 12
    stopping: string, optional,
            Sets the method used to stop the sifting process.
            None: Standard EMD equation .....
            TBA1: Second standard EMD equation ....
            resoultion: comes from ref [3]_. Need to set the parmeter res!
    shifiting_distance : float, optional
            Sets the minimum variance between IMF iterations.
            Default : 0.2
    res : float, optional
        stuff from ref [3]_      it is in dB
    Returns
    -------
    IMFs : ndarray
            An array of shape (len(data),N) where N is the number of found IMFs

    Notes
    -----

    References
    ----------
    .. [1] Huang H. et al. 1998 'The empirical mode decomposition and the
    Hilbert spectrum for nonlinear and non-stationary time series analysis.'
    Procedings of the Royal Society 454, 903-995

    .. [2] Zhao J., Huang D. 2001 'Mirror extending and circular spline
    function for empirical mode decomposition method'
    Journal of Zhejiang University (Science) V.2, No.3,P247-252

    .. [3] Rato R.T., Ortigueira M.D., Batista A.G 2008 'On the HHT,
    its problems, and some solutions.'
    Mechanical Systems and Signal Processing 22 1374-1394
    """

    #Set up signals array and IMFs array based on type of extrapolation
    # No extrapolation and 'extend' use signals array which is len(data)
    # Mirror extrapolation (Zhao 2001) uses a signal array len(2*data)
    if extrapolation == 'mirror':
        #Set up base
        base = len(data)
        nimfs = range(nimfs) # Max number of IMFs
        IMFs = np.zeros([base, len(nimfs)])
        ncomp = 0
        residual = data
        #Signals is 2*base
        signals = np.zeros([base*2, 2])
        #Mirror Dataset
        signals[0:base / 2, 0] = data[::-1][base / 2:]
        signals[base / 2:base + base / 2, 0] = data
        signals[base + base / 2:base * 2, 0] = data[::-1][0:base / 2]
        # Redfine base as len(signals) for IMFs
        base = len(signals)
        data_length = len(data) # Data length is used in recovering input data
        #Do spline fitting with periodic bounds
        periodic = 1
    else:
        base = len(data)
        signals = np.zeros([base, 2])
        nimfs = range(nimfs) # Max number of IMFs
        IMFs = np.zeros([base, len(nimfs)])
        ncomp = 0
        residual = data
        signals[:, 0] = data
        #Don't do spline fitting with periodic bounds
        periodic = 0

    for j in nimfs:
        # Extract at most nimfs IMFs no more IMFs to be found if Finish is True
        k = 0
        sd = 10
        finish = False

        while sd > shifting_distance and not(finish):

            #EMD magic here.
            min_arr, max_arr = exterma(signals[:,0])

            if min_arr.shape[1] <= 2 or max_arr.shape[1] <= 2:
                #If this IMF has become a straight line
                finish = True
            else:
                mean = envelope(min_arr, max_arr, base, n, periodic)

                if not(pearson):
                    alpha = 1
                else:
                    alpha = pearsonr(signals[:,0],mean)[0]
                signals[:,1] = signals[:,0] - alpha*mean

                #Calculate the shifting distance which is a measure of
                #simulartity to previous IMF
                if k > 0:
                        sd = np.sum((np.abs(signals[:,0] - signals[:,1])**2)) / np.sum(signals[:,0]**2)

                signals = signals[:,::-1]
                k += 1

            if finish:
                #If IMF is a straight line we are done here.
                IMFs[:,j]= residual
                ncomp += 1
                break

        if extrapolation == 'mirror':
            IMFs[:,j] = signals[data_length / 2:data_length
                                                           + data_length / 2,0]
            residual = residual - IMFs[:,j]#For j==0 residual is initially data

            #Mirror case requires IMF subtraction from data range then
            # re-mirroring for each IMF
            signals[0:data_length / 2,0] = residual[::-1][data_length / 2:]
            signals[data_length / 2:data_length + data_length / 2,0] = residual
            signals[data_length
                + data_length / 2:,0] = residual[::-1][0:data_length / 2]
            ncomp += 1

        else:
            IMFs[:,j] = signals[:,0]
            residual = residual - IMFs[:,j]#For j==0 residual is initially data
            signals[:,0] = residual
            ncomp += 1

    return IMFs