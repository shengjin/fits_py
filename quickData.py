# coding: utf-8
from astropy.io import fits
data = fits.getdata('LkCa15.continuum_map_r0.5.image.fits')
data = fits.getheader('LkCa15.continuum_map_r0.5.image.fits')
data = fits.getdata('LkCa15.continuum_map_r0.5.image.fits', ext=0)
data = fits.getdata('LkCa15.continuum_map_r0.5.image.fits', extname='primary')
