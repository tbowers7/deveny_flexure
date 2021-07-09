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


def extract_spectrum(spectrum,traces,nspix):
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
    """

    # Silence RuntimeWarning for overflow, this function only
    warnings.simplefilter('ignore', RuntimeWarning)

    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3
    return y


def find_lines(image, thresh=20.):
    """Automatically find and centroid lines in a 1-row image
 
    :image:
    :thresh: 20 DN above background

    Returns:
    :centers: List of line centers (pixel #)
    :fwhm: The computed FWHM
    """
    # Silence OptimizeWarning, this function only
    warnings.simplefilter('ignore', optimize.OptimizeWarning)

    _, ny = image.shape
    avgj = np.ndarray.flatten(image)

    # Create empty lists to fill
    peaks = []
    fwhm = []
 
    # Create background from median value of the image:
    bkgd = np.median(image)
    print(f'  Background level: {bkgd:.1f}')

    # Step through the cut and identify peaks:
    fmax = 11     # Minimum separations
    fwin = 15
    fhalfwin = int(np.floor(fwin/2))
    findmax = 50  # Maximum # of lines to find
    j0 = 0

    # Loop through lines
    for j in range(ny):
        
        if j > (ny - fmax):
            continue

        if avgj[j] > (bkgd + thresh):
            j1 = j
            if np.abs(j1 - j0) < fmax:      
                continue

            for jf in range(findmax):
                itmp0 = avgj[jf + j]
                itmp1 = avgj[jf + j + 1]
                if itmp1 < itmp0:
                    icntr = jf + j
                    break

            if (icntr < fmax/2) or (icntr > (ny - fmax/2 - 1)):
                continue

            xx = np.arange(fwin, dtype=float) + float(icntr - fhalfwin)
            temp = avgj[(icntr) - fhalfwin : (icntr) + fhalfwin + 1]
            temp = signal.medfilt(temp, kernel_size=3)
            
            # Run the fit, with error checking
            try:
                p0 = [1000, np.mean(xx), 3, bkgd]
                aa, _ = optimize.curve_fit(gaussfit_func, xx, temp, p0=p0)
            except RuntimeError:
                continue  # Just skip this one

            _ = gaussfit_func(xx, *aa)
            center = aa[1]
            fw = aa[2] * 2.355   # sigma -> FWHM
            if fw > 1.0:
                peaks.append(center)
                fwhm.append(fw)
            j0 = jf + j

    centers = np.asarray(peaks)
    
    cc = np.where(np.logical_and(centers >=0, centers <=2100))
    centers = centers[cc]

    print(f" Number of lines: {len(centers)}")
 
    return (centers, fwhm)


def specavg(spectrum, trace, wsize):
    """Extract an average spectrum along trace of size wsize
    :param spectrum: input spectrum
    :param trace: the trace along which to extract
    :param wsize: the size of the extraction (usually odd)
    :return:
    """
    # Case out the dimensionality of traces... 0 -> return
    if spectrum.ndim == 0:
        return 0
    elif spectrum.ndim == 1:
        nx = spectrum.size
    else:
        _, nx = spectrum.shape
    speca = np.empty(nx, dtype=float)

    whalfsize = int(np.floor(wsize/2))

    # Because of python indexing, we need to "+1" the upper limit in order
    #   to get the full wsize elements for the average
    for i in range(nx):
        speca[i] = np.average(spectrum[int(trace[i]) - whalfsize : 
                                       int(trace[i]) + whalfsize + 1, i])
    
    return speca.reshape((1,nx))
