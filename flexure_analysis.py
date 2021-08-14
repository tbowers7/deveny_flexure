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
    """flexure_analysis Driving routine for the analysis

    [extended_summary]

    Parameters
    ----------
    data_dir : `str`
        Directory where the data live
    """    

    for grating in ['DV1','DV2','DV5']:

        # Create an ImageFileCollection with files matching this grating
        gcl = load_images(data_dir, grating)

        # If no files for this grating, move along
        if len(gcl.files) == 0:
            continue

        # AstroPy Table of line positions for each image in the IFC
        table = get_line_positions(gcl)
        # Go through the identified lines and produce a set found in all images
        table = validate_lines(table)
        table.pprint()

        # Write out the table to disk
        with open(f"line_positions_{grating}.csv",'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(table)

        make_plots(table, grating)

    return


def load_images(data_dir, grating):
    """load_images Load in the images associated with DATA_DIR and grating
    
    [extended_summary]

    Parameters
    ----------
    data_dir : `str`
        The directory containing the data to analyze
    grating : `str`
        The grating ID to use

    Returns
    -------
    `ccdproc.image_collection.ImageFileCollection`
        IFC of the files meeting the input criteria
    """    
    # Dictionary
    gratid = {'DV1':'150/5000', 'DV2':'300/4000', 'DV5':'500/5500'}

    # Load the images of interest in to an ImageFileCollection()
    icl = ccd.ImageFileCollection(data_dir)

    # Return an ImageFileCollection filtered by the grating desired
    return icl.filter(grating=gratid[grating])


def get_line_positions(icl, win=11, thresh=10000.):
    """get_line_positions Compute the line positions for the images in the icl

    [extended_summary]

    Parameters
    ----------
    icl : `ccdproc.image_collection.ImageFileCollection`
        ImageFileCollection of images to work with
    win : `int`, optional
        Window (in pixels) across which to extract the spectrum, [Default: 11]
    thresh : `float`, optional
        Line intensity (ADU) threshold for detection, [Default: 10000.]

    Returns
    -------
    `astropy.table.table.Table`
        Table of line positions with associated metadata
    """    
    # Put everything into a list of dicionaties
    flex_line_positions = []

    # This will only give the x values of the fits file.
    # For each of the images,
    for ccd, fname in icl.ccds(return_fname=True):
        # Need a lower threshold for DV5 than for DV1
        if ccd.header['grating'] == '500/5500':
            thresh = 1000.
        # Check for bias frames
        if ccd.header['exptime'] == 0:
            continue
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
        print(f"Found {nc} Line Centers: {[f'{cent:.1f}' for cent in centers]}")
        #====================

        h = ccd.header
 
        flex_line_positions.append({'filename':fname,
                                    'obserno': h['obserno'],
                                    'alt':h['telalt'],
                                    'az':h['telaz'],
                                    'cass':h['rotangle'],
                                    'utcstart':h['utcstart'],
                                    'lamps':h['lampcal'],
                                    'grating':h['grating'],
                                    'grangle':h['grangle'],
                                    'slitwidth':h['slitasec'],
                                    'nlines':nc,
                                    'xpos':centers})

    t = Table(flex_line_positions)
    return t


def validate_lines(t):
    """validate_lines Validate the found lines to produce a uniform set

    The number of lines identified will vary form image to image.  This
    function validates the lines to return the set of lines found in ALL
    images for this grating.

    Parameters
    ----------
    t : `astropy.table.table.Table`
        AstroPy Table as produced by get_line_positions()

    Returns
    -------
    `astropy.table.table.Table`
        AstroPy Table identical to input except the lines are validated
    """    
    return t


def make_plots(t, grating):
    """make_plots Make plots of the data... with subcalls to fitting functions

    [extended_summary]

    Parameters
    ----------
    t : `astropy.table.table.Table`
        AstroPy Table for this grating, containing all the data!
    grating : `str`
        Grating name, for labeling plots and creating filenames
    """    

    # Set up the plotting environment
    _, ax = plt.subplots()
    tsz = 8

    

    ax.set_xlabel('Something', fontsize=tsz)
    ax.set_ylabel('Something', fontsize=tsz)

    # Final adjustments and save figure
    ax.tick_params('both', labelsize=tsz, direction='in', top=True, right=True)
    plt.tight_layout()
    plt.savefig(f"flexure_analysis_{grating}.eps")
    plt.savefig(f"flexure_analysis_{grating}.png")

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
