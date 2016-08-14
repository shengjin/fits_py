# coding: utf-8
from astropy.io import fits
hdulist = fits.open('LkCa15.continuum_map_r0.5.image.fits')
hdulist.info()
hdu = hdulist[0]
hdu.data.shape
hdu.header
hdu.header['TELESCOP']
hdu.header['INSTRUME']
plt.imshow(hdu.data[0,0,:,:], origin='lower')
import matplotlib.pyplot as plt
plt.imshow(hdu.data[0,0,:,:], origin='lower')
plt.show()
hdu.header.remove('CRPIX3')
hdu.header.remove('CRVAL3')
hdu.header.remove('CDELT3')
hdu.header.remove('CUNIT3')
hdu.header.remove('CTYPE3')
hdu.header.remove('CTYPE4')
hdu.header.remove('CUNIT4')
hdu.header.remove('CDELT4')
hdu.header.remove('CRVAL4')
hdu.header.remove('CRPIX4')
hdu.writeto('lat_background_model_slice.fits')
hdulist.writeto('lat_background_model_slice_allhdus.fits')
