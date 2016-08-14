#!/usr/bin/env python2
from astropy.io import fits
import matplotlib.pyplot as plt
hdulist = fits.open('LkCa15_12CO_cube.image.fits')
hdulist.info()
hdu = hdulist[0]
m,n,o,p=hdu.data.shape
hdu.header
print m,n,o,p
#plt.imshow(hdu.data[0,0,:,:], origin='lower')
#plt.show()
#hdu.writeto('lat_background_model_slice.fits')
#hdulist.writeto('lat_background_model_slice_allhdus.fits')




ImageJyppix_scaled_ave=hdu.data[0,0,:,:]
for i in range(n-1):
    ImageJyppix_scaled_ave= ImageJyppix_scaled_ave+hdu.data[0,i+1,:,:]

plt.imshow(ImageJyppix_scaled_ave[:,:], origin='lower')
#plt.xlabel('AU')
#plt.ylabel('AU')
#plt.ylim(-180,180)
#plt.xlim(-180,180)
#cbar1=plt.colorbar()
#cbar1.set_label("Janskey/pixel")
plt.title("Flux density")
plt.savefig('flux_density_ave.png')
plt.clf()

