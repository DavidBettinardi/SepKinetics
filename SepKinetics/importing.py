#!/bin/env python2.7

import sys
import pandas as pd
import matplotlib.pyplot as plt
import plotly.offline as py
import plotly.graph_objs as go
import peakutils

'''# TODO:
-check Martin's disseration for deconvolution
-add more tests
-add google style comments to code
-import .xlsx instead of .csv (comes from clarity as .xlsx)?
-save 3D plots?
-fix peakfinding x-axis
-add 2D peak labels
-fix 3D y-axis
-add plotting series of histograms (see PythonProject folder in chrome)
-create PyPi account for distrubtion
'''

""" 04-22-18 - David j Bettinardi - SepKinetics software
The script imports .csv kinetic data from the Olis CLARiTY
spectrophotometer and outputs 2D and 3D visualizers.
Three datastructure approaches were attempted; current approach
uses dataframes in -m pandas"""

'''call in terminal using $ python importing.py filename.csv'''

#filename = 'kintest.csv'    # for test_functions, overridden in main

# pandas approach [third try]
def main(file):
    """Creates dataframe with top row as headers"""
    try:
        return df
    except:
        df = pd.read_csv(file, header = 0) # first row as header
    return df


def two_dplot(specdata):
    """Makes an overlapping 2D plot of wavelength vs. kinetics @ all timepoints"""
    if not isinstance(specdata, pd.DataFrame):
        specdata = main(specdata)
        specdata.set_index('0', inplace = True)   #reset index as time-series
    else:
        pass
    return specdata.plot(legend=False, colormap = 'plasma')


def three_dplot(specdata):
    """Makes a 3D surface plot"""
    if not isinstance(specdata, pd.DataFrame):
        specdata = main(specdata)
        specdata.set_index('0', inplace = True)   #reset index as time-series
    else:
        pass
    #specdata2 = specdata.copy()
    #specdata2.set_index('0', inplace = True)

    #y = specdata.iloc[:,0].tolist()   # attempt to overwrite y-axis

    '''Working plotter, wrong y-axis'''    # MUST FIX
    data = [go.Surface(z=specdata.as_matrix(), colorscale='Viridis')]     # transform to matrix for surface plot

    layout = go.Layout(
                    scene = dict(
                    xaxis = dict(
                        title='Time (s)'),
                    yaxis = dict(
                        title='Wavelength (nm)'),
                    zaxis = dict(
                        title='Abs'),),
                    width=700,
                  )
    fig = go.Figure(data=data, layout=layout)
    return py.plot(fig)


if __name__ == '__main__':
    try:
        specdata = sys.argv[1]
    except:
        raise ValueError('Specify the .csv you wish to import :)')

    #two_dplot(specdata)
    #three_dplot(specdata)
    #plt.show()
