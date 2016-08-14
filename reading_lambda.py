#!/usr/bin/env python2
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

plot = False

cc = 29979245800.0       # Speed of light (cm/s)

#hdulist = fits.open('LkCa15.continuum_map_r0.5.image.fits')
#hdulist = fits.open('HLTau_B7cont_av20_mscale_ap_shifted.image.fits')
hdulist = fits.open('LkCa15_12CO_cube.image.fits')

hdulist.info()
hdu = hdulist[0]


naxis3=hdu.header['NAXIS3']
crval3=hdu.header['CRVAL3']
cdelt3=hdu.header['CDELT3']
crpix3=hdu.header['CRPIX3']

grid = np.arange(naxis3) + 1.0
freq = crval3 + (grid-crpix3) * cdelt3
lamb = cc/freq*1e4

np.savetxt('lamb_alma', np.transpose([lamb]))

#print hdu.header
if plot:
    plt.imshow(hdu.data[0,0,:,:], origin='lower')
    plt.colorbar()
    plt.ylim(300,724)
    plt.xlim(300,724)
    plt.savefig('dust.png')
    plt.clf()

