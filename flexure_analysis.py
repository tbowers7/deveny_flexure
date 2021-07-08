import csv
import ccdproc as ccd
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt


def flexure_analysis(args):
    # Driving routiune for the analysis

    for grating in ['DV1','DV2','DV5']:
        gcl = load_images(grating)
        summary = gcl.summary['obserno','telalt','telaz','rotangle'].pprint()
        table = get_line_positions(gcl)
        table.pprint()
        with open('table_containing_line_positions_DV1.csv','w', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(table)

    make_plots()
    return


def load_images(grating):
    # Load in the Images of interest

    icl = ccd.ImageFileCollection('./data/all_data')
    DV1 = '150/5000'
    DV2 = '300/4000'
    DV5 = '500/5500'
#    gid = '150/5000' if grating == 'DV1' else '500/5500' # do I need to do anything for the other gratings?
    return icl.filter(grating=DV1)

def find_indices(arr,condition):
    return [i for i, elem in enumerate(arr) if condition(elem.all())]

def get_line_positions(icl):
    """
    Inputs: 'icl' = ImageFileCollection of images to work with
    """
    # Put everything into a list of dicionaties
    flex_line_positions = []

    # This will only give the x values of the fits file.
    # For each of the images,
    for ccd, fname in icl.ccds(return_fname=True):
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
                                 'Az':somefunctionthatcalssto('telaz'),
                                 'rot':somefunctionthatcallsto('rotangle)})
            # save lines to cached lines array
            # for each scan: if above flat threshold,
                # mark center value of line,
                # find start of increase above flat level
                # find end of increase above flat level
                # define limits of line
                # find middle of the line
                # count that x value as a valid line and store
        print(lines_by_row[-1])
        flex_line_positions.append({'file':fname,
                                    'lines': lines_by_row})
    t = Table(flex_line_positions)
    return t



def make_plots():
    # Make plots of the data... with subcalls to fitting functions
    return





def main(args):
    flexure_analysis(args)
    return

if __name__ == '__main__':
    import sys
    main(sys.argv)
