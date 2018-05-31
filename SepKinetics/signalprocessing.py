#!/bin/env python2.7

from importing import main, two_dplot, three_dplot
import sys
import pandas as pd
import matplotlib.pyplot as plt
import plotly.offline as py
import plotly.graph_objs as go
import peakutils
import numpy as np
#from sklearn.decomposition import PCA


def baseline(specdata):
    """For a single spectrum: Performs a baseline-identifying algorithm
    from -m peakutils, corrects, and plots"""
    #indices, midpoint_spectrum = peakfind() # obtain midpoint peak data

    # execute the baseline-identifier on this spectrum
    baseline_values = peakutils.baseline(specdata)

    # print the correction for validation
    original_spectrum = go.Scatter(      # show original spectrum
        x=[j for j in range(len(specdata))],
        y=specdata,
        mode='lines',
        name='Abs spectrum of the midpoint-time in kinetic series'
    )

    baseline = go.Scatter(      # show baseline
        y=baseline_values,
        x=[j for j in range(len(specdata))],
        mode='markers',
        marker=dict(
            size=3,
            color='rgb(255,0,0)',
            symbol='ope-circle'
        ),
        name='Basline values'
    )

    corrected_spectrum = go.Scatter(        # show baseline correction
        y=specdata - baseline_values,
        x=[j for j in range(len(specdata))],
        mode='lines',
        marker=dict(
            color='rgb(155,100,0)',
        ),
        name='Basline subtracted'
    )

    data = [original_spectrum, baseline, corrected_spectrum]     # overlay all three
    return py.plot(data)


def series_baseline(seriesdata):
    """For a spectral series: Performs a baseline-identifying algorithm
    from -m peakutils, and corrects"""

    df = main(seriesdata)
    df.set_index('0', inplace = True)
    baseline_series = df.copy()     # make a deep copy to the dataframe
    corrected_spectra = df.copy()   # make a deep copy to the dataframe

    for i in xrange(len(df.iloc[0,:])):     # execute baseline-identifier algorithm by column
        baseline_series.iloc[:,i] = peakutils.baseline(df.iloc[:,i], deg = 3)    # overwrite the df data
    corrected_spectra = df - baseline_series    # subtract from original spectra to find correction
    corrected_spectra_indexes = corrected_spectra.reset_index()
    #print corrected_spectra_indexes
    return baseline_series, corrected_spectra, corrected_spectra_indexes


def series_peakfind(specdata, threshold = 0.2, min_dist=50):
    "Produces a list of lists containing the peaks at each timepoint"
    if not isinstance(specdata, pd.DataFrame):  # check datatype
        specdata = main(specdata)
        specdata.set_index('0', inplace = True)   #reset index as time-series
        print 'Warning: This series is not baseline corrected!'
    else:
        pass

    indexes=[]      # initialize the empty array
    for i in xrange(len(specdata.iloc[0,:])):     # execute peakfinding algorithm by column
        individual_spectrum = specdata.iloc[:,i]
        indexes.append(peakutils.indexes(individual_spectrum, min_dist= min_dist, thres = threshold/max(individual_spectrum)))
    return indexes


def midpoint(specdata):
    """Finds major/minor peaks according to threshold values; Threshold
    value is divided by max peak height"""
    if not isinstance(specdata, pd.DataFrame):  # check datatype
        specdata = main(specdata)
        specdata.set_index('0', inplace = True)   #reset index as time-series
    else:
        pass
    midpoint = (len(specdata.iloc[0,:]))/2    # find the mid-timepoint
    midpoint_spectrum = specdata.iloc[:, midpoint]    # find the midpoint spectrum
    return midpoint_spectrum


def peakdet(v, delta, x = None):
    """
    From: https://gist.github.com/endolith/250860
    """
    maxtab = []
    mintab = []

    if x is None:
        x = np.arange(len(v))

    v = np.asarray(v)
    print v

    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')

    if not np.isscalar(delta):
        sys.exit('Input argument delta must be a scalar')

    if delta <= 0:
        sys.exit('Input argument delta must be positive')

    mn, mx = np.Inf, -np.Inf
    mnpos, mxpos = np.NaN, np.NaN

    lookformax = True

    for i in np.arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return np.array(maxtab)


def midpoint_peakfind(specdata):
    """Finds major peaks according to threshold values; Threshold
    value is divided by max peak height"""

    if not isinstance(specdata, pd.DataFrame):  # check datatype
        specdata = main(specdata)
        specdata.set_index('0', inplace = True)   #reset index as time-series
    else:
        pass

    midpoint = (len(specdata.iloc[0,:]))/2    # find the mid-timepoint
    midpoint_spectrum = specdata.iloc[:, midpoint]    # find the midpoint spectrum

    # exectue the peak finding on this spectrum
    indices = signal.find_peaks_cwt(midpoint_spectrum)
        # trace the midpoint spectrum
    spectrum = go.Scatter(
        x=[j for j in range(len(midpoint_spectrum))],
        y=midpoint_spectrum,
        mode='lines',
        name='Abs spectrum of the midpoint-time in kinetic series'
    )

    # mark the peaks
    peaks = go.Scatter(
        x=indices,
        y=[midpoint_spectrum.iloc[j] for j in indices],
        mode='markers',
        marker=dict(
            size=8,
            color='rgb(255,0,0)',
            symbol='cross'
        ),
        name='Detected Peaks'
    )

    data = [spectrum, peaks]
    return py.plot(data)


def initial_rates(corrected, time, lam1=88, lam2=110, lam3=125):
    """Takes a kinetic dataset and the initial rate method for
    determining kinetics constants for zero-order reactions at
    over a period of time at the specified wavelengths (default
    = Nd peaks)"""
    print 'Calculating Rate Constants by Initial Rate Method...'

    lambda1=(corrected.iloc[lam1,15:])
    lambda2=(corrected.iloc[lam2,15:])
    lambda3=(corrected.iloc[lam3,15:])

    t1 = [float(i) for i in lambda1.index.values.tolist()]  # time is shared among spectra

    regress1 = np.polyfit(t1[0:time], lambda1[0:time], 2)    # find a polynomial fit
    regress2 = np.polyfit(t1[0:time], lambda2[0:time], 2)    # find a polynomial fit
    regress3 = np.polyfit(t1[0:time], lambda3[0:time], 2)    # find a polynomial fit

    p1 = np.poly1d(regress1)  # extablish the standalone polynomial eq
    p2 = np.poly1d(regress2)  # extablish the standalone polynomial eq
    p3 = np.poly1d(regress3)  # extablish the standalone polynomial eq

    plist1 = [p1(i) for i in t1[0:time]]  # find value at each timepoint
    plist2 = [p2(i) for i in t1[0:time]]  # find value at each timepoint
    plist3 = [p3(i) for i in t1[0:time]]  # find value at each timepoint

    # plot data and regressors together
    plt.plot(t1[0:(time*2)], lambda1[0:(time*2)])
    plt.plot(t1[0:(time*2)], lambda2[0:(time*2)])
    plt.plot(t1[0:(time*2)], lambda3[0:(time*2)])
    plt.plot(t1[0:time], plist1, '--')
    plt.plot(t1[0:time], plist2, '--')
    plt.plot(t1[0:time], plist3, '--')
    plt.legend(title='Wavelengths')
    plt.title('Absorption kinetics with overlayed regressors')
    plt.show()

    d1 = np.polyder(p1)   # find the derivative function of the polynomial
    d2 = np.polyder(p2)   # find the derivative function of the polynomial
    d3 = np.polyder(p3)   # find the derivative function of the polynomial

    # print the derivative of the regressor at t=0 to find the initial rate
    print 'k_observed @', corrected.index[lam1], 'nm = ', d1(t1[0])
    print 'k_observed @', corrected.index[lam2], 'nm = ', d2(t1[0])
    print 'k_observed @', corrected.index[lam3], 'nm = ', d3(t1[0])

    return d1(t1[0]), d2(t1[0]), d3(t1[0])


if __name__ == '__main__':
    try:

        specdata = sys.argv[1]

    except:
        raise ValueError("Specify the .csv you wish to import")

    baseline_series, corrected_spectra, corrected_spectra_indexes = series_baseline(specdata)
    #two_dplot(specdata), two_dplot(baseline_series), two_dplot(corrected_spectra)
    #two_dplot(corrected_spectra.iloc[50:150,1:])
    #plt.show()
    #three_dplot(corrected_spectra)

    initial_rates(corrected_spectra, 100)

    '''
    k1a= (lambda1[20] - lambda1[0])/(t1[20]-t1[0])
    k1b= (lambda1[40] - lambda1[20])/(t1[20]-t1[0])
    k1c= (lambda1[60] - lambda1[40])/(t1[20]-t1[0])
    if k1a<k1b or k1b<k1c:
        print 'bad slopes'

    print k1a, k1b, k1c

    difference = ((k1b-k1a) + (k1c-k1b))/2
    print 'k(0) @576nm = ', k1a - difference

    k2a= (lambda2[20] - lambda2[0])/(t1[20]-t1[0])
    k2b= (lambda2[40] - lambda2[20])/(t1[20]-t1[0])
    k2c= (lambda2[60] - lambda2[40])/(t1[20]-t1[0])
    difference = ((k1b-k1a) + (k1c-k1b))/2
    print 'k(0) @583nm = ', k2a - difference

    k3a= (lambda3[20] - lambda3[0])/(t1[20]-t1[0])
    k3b= (lambda3[40] - lambda3[20])/(t1[20]-t1[0])
    k3c= (lambda3[60] - lambda3[40])/(t1[20]-t1[0])
    difference = ((k1b-k1a) + (k1c-k1b))/2
    print 'k(0) @588nm = ', k3a - difference
    '''


    ''' Peak finding util
    from matplotlib.pyplot import plot, scatter, show
    baseline_series, corrected_spectra = series_baseline(specdata)
    midpoint = corrected_spectra.iloc[60:150,199]
    series = np.asarray(midpoint)
    maxtab = peakdet(series,.0015)
    plot(series)
    scatter(np.array(maxtab)[:,0], np.array(maxtab)[:,1], color='red')
    show()
    '''




    # Implementation of Principal Component Analysis (PCA)
    #specdata = main(specdata)
    #specdata.set_index('0', inplace = True)
    #print specdata
    #array = specdata.values
    #print array
    #pca = PCA(n_components = 1)
    #pca.fit(array)
    #comp = pca.components_
    #print comp
    #pca = PCA().fit(array)
    #plt.plot(np.cumsum(pca.explained_variance_ratio_))
    #plt.xlabel('number of components')
    #plt.ylabel('cumulative explained variance')
    #plt.show()
