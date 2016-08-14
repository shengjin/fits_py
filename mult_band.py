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
#plt.imshow(hdu.data[0,i*2,:,:], origin='lower')
#plt.show()
#hdu.writeto('lat_background_model_slice.fits')
#hdulist.writeto('lat_background_model_slice_allhdus.fits')


n_row=5
n_col=9

plt, axes = plt.subplots(nrows=5, ncols=9, figsize=(13, 8), dpi=80, sharex=True, sharey=True)


for i in range(n_row):
    for j in range(n_col):
        print i,j
        im=axes[i,j].imshow(hdu.data[0,(n_col*i+j)+22,:,:], origin='lower')
        axes[i,j].set_ylim(300,724)
        axes[i,j].set_xlim(300,724)
        #im=axes[i,j].pcolormesh(im_x_au, im_y_au, ImageJyppix_scaled[:,:,n_col*i+j], vmin=v_min, vmax=v_max, cmap='jet')
        #im=axes[i,j].pcolormesh(im_x_au, im_y_au, ImageJyppix_scaled[:,:,n_col*i+j], vmin=v_min, vmax=v_max, cmap='RdBu')
        #im=axes[i,j].pcolormesh(im_x_au, im_y_au, (ImageJyppix_scaled[:,:,n_col*i+j]-Image_dust_Jyppix_scaled[:,:,0]), vmin=v_min, vmax=v_max, cmap='RdBu')


plt.tight_layout()

plt.subplots_adjust(hspace=0.0, wspace=0.0)
plt.savefig('gas_lines_radmc.png')

plt.clf()











