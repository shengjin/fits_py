# coding: utf-8
from astropy.io import fits
import numpy as np
hdu = fits.PrimaryHDU()
hdu.data = np.random.random((128,128))
hdu.header['telescop'] = 'Python Observatory'
hdu.writeto('random_array.fits')
hdu.writeto('random_array.fits', clobber=True) # overwrite
