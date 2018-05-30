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
    """Finds major/minor peaks according to threshold values; Threshold
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

    lambda1=-np.log(corrected_spectra.iloc[88,15:])
    lambda2=-np.log(corrected_spectra.iloc[110,15:])
    lambda3=-np.log(corrected_spectra.iloc[125,15:])

    lambda1.plot(legend=True)
    lambda2.plot(legend=True)
    lambda3.plot(legend=True)

    #plt.show()

    x = [float(i) for i in lambda1.index.values.tolist()]
    regress = np.polyfit(lambda1,x,3)
    p = np.poly1d(regress)
    d = np.polyder(p)
    print p
    print d




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
