#!/bin/env python2.7

"""Test scipts for importing.py"""

from .. import importing
import numpy as np
import os

def test_main():
    fn = os.path.join(os.path.dirname(__file__), 'kintest.csv')
    df = importing.main(fn)
    df.set_index('0', inplace = True)
    timepoints = [float(i) for i in list(df)]

    assert len(timepoints) == 2
    np.testing.assert_almost_equal(df.iloc[0,0], 0.004224, decimal=5)

def test_three_dplot():
    pass
