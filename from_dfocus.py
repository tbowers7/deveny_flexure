# -*- coding: utf-8 -*-
#
#  This file is part of ______.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 08-Jul-2021
#
#  @author: tbowers, bshafransky

"""Selected (trimmed) routines from `dfocus`

These are from the pydeveny.dfocus() code, pruned for the immediate use case.
"""

import numpy as np
import warnings
from scipy import optimize
from scipy import signal


def extract_spectrum(spectrum, traces, nspix):
    """Object spectral extraction routine
   
    Input:
      spectrum: 2D spectral image
      traces: Trace line(s) along which to extract the spectrum
      nspix: Window width across which to extract the spectrum

    Output:
      spectra: 2-d array of spectra of individual orders
    """
    
    # Set # orders, size of each order based on traces dimensionality; 0 -> return
    if traces.ndim == 0:
        return 0                
    norders, nx = (1, traces.size) if traces.ndim == 1 else traces.shape

    # Start out with an empty array
    spectra = np.empty((norders, nx), dtype=float)

    # Get the averaged spectra
    for io in range(norders):
        spectra[io,:] = specavg(spectrum, traces[io,:], nspix)

    return spectra


def gaussfit_func(x, a0, a1, a2, a3):
    """Simple Gaussian function for fitting line profiles
    Input:
      x: array of x values
      a0: Amplitude
      a1: Mean x value (mu)
      a2: Gaussian sigma
      a3: Background offset

    Output:
      y: array of y values corresponding to input a's
    """

    # Silence RuntimeWarning for overflow, this function only
    warnings.simplefilter('ignore', RuntimeWarning)

    z = (x - a1) / a2
    return a0 * np.exp(-z**2 / 2.) + a3


def find_lines(image, thresh=20., findmax=50, minsep=11, fit_window=15):
    """Automatically find and centroid lines in a 1-row image
 
    :image:
    :thresh: Threshold above which to indentify lines [Default: 20 DN above bkgd]
    :findmax: Maximum number of lines to find [Default: 50]
    :minsep: Minimum line separation for identification [Default: 11 pixels]
    :fit_window: Size of the window to fit Gaussian [Default: 15 pixels]

    Returns:
    :centers: List of line centers (pixel #)
    :fwhm: The computed FWHM
    """
    # Silence OptimizeWarning, this function only
    warnings.simplefilter('ignore', optimize.OptimizeWarning)

    # Define the half-window
    fhalfwin = int(np.floor(fit_window/2))
 
     # Get size and flatten to 1D
    _, nx = image.shape
    spec = np.ndarray.flatten(image)

    # Find background from median value of the image:
    bkgd = np.median(spec)
    print(f'  Background level: {bkgd:.1f}')

    # Create empty lists to fill
    cent, fwhm = ([], [])
    j0 = 0

    # Step through the cut and identify peaks:
    for j in range(nx):
        
        # If we get too close to the end, skip
        if j > (nx - minsep):
            continue

        # If the spectrum at this pixel is above the THRESH...
        if spec[j] > (bkgd + thresh):

            # Mark this pixel as j1
            j1 = j

            # If this is too close to the last one, skip
            if np.abs(j1 - j0) < minsep:      
                continue

            # Loop through 0-FINDMAX...  (find central pixel?)
            for jf in range(findmax):
                itmp0 = spec[jf + j]
                itmp1 = spec[jf + j + 1]
                if itmp1 < itmp0:
                    icntr = jf + j
                    break

            # If central pixel is too close to the edge, skip
            if (icntr < minsep/2) or (icntr > (nx - minsep/2 - 1)):
                continue

            # Set up the gaussian fitting for this line
            xmin, xmax = (icntr - fhalfwin, icntr + fhalfwin + 1)
            xx = np.arange(xmin, xmax, dtype=float)
            temp = spec[xmin : xmax]
            # Filter the SPEC to smooth it a bit for fitting
            temp = signal.medfilt(temp, kernel_size=3)
            
            # Run the fit, with error checking
            try:
                p0 = [1000, np.mean(xx), 3, bkgd]
                aa, _ = optimize.curve_fit(gaussfit_func, xx, temp, p0=p0)
            except RuntimeError:
                continue  # Just skip this one

            # If the width makes sense, save
            if (fw := aa[2] * 2.355) > 1.0:      # sigma -> FWHM
                cent.append(aa[1])
                fwhm.append(fw)

            # Set j0 to this pixel before looping on
            j0 = jf + j

    # Make list into an array, check again that the centers make sense
    centers = np.asarray(cent)    
    c_idx = np.where(np.logical_and(centers > 0, centers <= nx))
    centers = centers[c_idx]

    print(f" Number of lines: {len(centers)}")
 
    return (centers, fwhm)


def specavg(spectrum, trace, wsize):
    """Extract an average spectrum along trace of size wsize
    :param spectrum: input spectrum
    :param trace: the trace along which to extract
    :param wsize: the size of the extraction (usually odd)
    :return:
    """
    # If ndim = 0, return, otherwise get nx
    if spectrum.ndim == 0:
        return 0
    nx = (spectrum.shape)[-1]

    speca = np.empty(nx, dtype=float)
    whalfsize = int(np.floor(wsize/2))

    # Because of python indexing, we need to "+1" the upper limit in order
    #   to get the full wsize elements for the average
    for i in range(nx):
        speca[i] = np.average(spectrum[int(trace[i]) - whalfsize : 
                                       int(trace[i]) + whalfsize + 1, i])
    
    return speca.reshape((1,nx))
