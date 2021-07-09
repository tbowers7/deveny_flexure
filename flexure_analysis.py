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

"""Analysis of the flexure seen in the DeVeny Spectrograph (LDT)

This file contains the main driver for the analysis.
Should be run in an environment containing:
    * AstroPy
    * CCDPROC
    * NumPy
    * Matplotlib

Run from the command line:
% python flexure_analysis.py DATA_DIR

"""

import csv
import ccdproc as ccd
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

# Get local routines
from from_dfocus import *


def flexure_analysis(data_dir):
    """Driving routine for the analysis
    
    Input:
        `data_dir` : Directory where the data are
    """

    for grating in ['DV1']: #['DV1','DV2','DV5']:

        # Create an ImageFileCollection with files matching this grating
        gcl = load_images(data_dir, grating)

        # Summary Table
        summary = gcl.summary['obserno','telalt','telaz','rotangle']

        # Data table of line positions for each image in the ImageFileCollection
        table = get_line_positions(gcl)
        table.pprint()

        # Write out the table to disk
        with open('table_containing_line_positions_DV1.csv','w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(table)

    make_plots()
    return


def load_images(data_dir, grating):
    """Load in the images associated with DATA_DIR and grating
    
    Inputs: `data_dir`: The directory containing the data to analyze
            `grating`: The grating ID to use
    """
    # Dictionary
    gratid = {'DV1':'150/5000', 'DV2':'300/4000', 'DV5':'500/5500'}

    # Load the images of interest in to an ImageFileCollection()
    icl = ccd.ImageFileCollection(data_dir)

    # Return an ImageFileCollection filtered by the grating desired
    return icl.filter(grating=gratid[grating])


def find_indices(arr,condition):
    return [i for i, elem in enumerate(arr) if condition(elem.all())]


def get_line_positions(icl, win=11, thresh=10000.):
    """
    Inputs: 'icl' = ImageFileCollection of images to work with

            These inputs are for the `dfocus`-derived code, and may be left
            at their default values.
            'win' = Something about the window to extract
            'thresh' = ADU threshold, above which look for lines
    """
    # Put everything into a list of dicionaties
    flex_line_positions = []

    # This will only give the x values of the fits file.
    # For each of the images,
    for ccd, fname in icl.ccds(return_fname=True):
        print("")
        #====================
        # Code cut-and-paste from dfocus() -- Get line centers above `thresh`
        # Parameters for DeVeny (2015 Deep-Depletion Device):
        n_spec_pix, prepix = (2048, 50)
        # Trim the image (remove top and bottom rows, pre- and post-scan pixels)
        spec2d = ccd.data[12:512,prepix:prepix+n_spec_pix]
        ny, nx = spec2d.shape
        trace = np.full(nx, ny/2, dtype=float).reshape((1,nx)) # Right down the middle
        spec1d = extract_spectrum(spec2d, trace, win)
        # Find the lines:
        centers, _ = find_lines(spec1d, thresh=thresh)
        nc = len(centers)
        print(F"In get_line_positions(), number of lines: {nc}")
        print(f"Line Centers: {[f'{cent:.1f}' for cent in centers]}")
        #====================


        #### Ben: At this point, you will have the line centers for this image,
        #         and you can then proceed to place those into a table for
        #         future analysis. 


        # detrmine flat level of image
        #take a defined set of known flat collums
        max_row0=max(np.transpose(ccd.data)[0])
        max_row1=max(np.transpose(ccd.data)[1])
        max_row2=max(np.transpose(ccd.data)[2])
        flat_value =  max(max_row0,max_row1,max_row2)

        # measure the line position
        # scan each row,
        lines_by_row = []

        for r in range(50,100):
            # if the majority of the collumns give a value above flat, store index in array
            lines_by_row.append({'row':r,
                                 'lines': find_indices(ccd.data[r] > flat_value, lambda e: e ),
                                 'alt':somefunctionthatcallsto('telalt'),
                                 'Az':somefunctionthatcallsto('telaz'),
                                 'rot':somefunctionthatcallsto('rotangle')})
            # save lines to cached lines array
            # for each scan: if above flat threshold,
                # mark center value of line,
                # find start of increase above flat level
                # find end of increase above flat level
                # define limits of line
                # find middle of the line
                # count that x value as a valid line and store
        #print(lines_by_row[-1])
        flex_line_positions.append({'file':fname,
                                    'lines': lines_by_row})
    t = Table(flex_line_positions)
    return t


def somefunctionthatcallsto(str):
    return 0


def make_plots():
    # Make plots of the data... with subcalls to fitting functions
    return




#==============================================================================
def main(args):
    """Main driving routine

    Call should be of form:
    % python flexure_analysis.py DATA_DIR
    """
    from os import path

    if len(args) == 1:
        print(f"ERROR: scrpit {args[0]} requires the DATA_DIR to analyze.")
        return
    if not path.isdir(args[1]):
        print(f"ERROR: DATA_DIR must be a directory containing the data to analyze.")
        return
    if len(args) > 2:
        print(f"WARNING: I'm ignoring the following arguments: {args[2:]}")

    # Run the analysis
    flexure_analysis(args[1])

    return


if __name__ == '__main__':
    import sys
    main(sys.argv)
