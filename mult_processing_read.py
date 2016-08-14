import os
import glob
import multiprocessing
import shutil

import pyfits
from scipy.ndimage import median_filter

# Define a function to run on files. The steps are:
# - read in FITS file
# - convolve the data in the primary HDU
# - write out the result to a new file
def smooth(filename):
    print "Processing %s" % filename
    hdulist = pyfits.open(filename)
    hdulist[0].data = median_filter(hdulist[0].data, 15)
    hdulist.writeto(filename.replace('files/', 'files_smooth/'),
                    clobber=True)

# Search for all FITS files
files = glob.glob('files/*.fits')

# Remove output directory if it already exists
if os.path.exists('files_smooth'):
    shutil.rmtree('files_smooth')

# Create output directory
os.mkdir('files_smooth')

# Define a 'pool' of 12 processes
p = multiprocessing.Pool(processes=12)

# Run the function over all files in parallel
result = p.map(smooth, files)
p = multiprocessing.Pool(processes=12)

%time result = p.map(smooth, files)

p = multiprocessing.Pool(processes=1)

%time result = p.map(smooth, files)


