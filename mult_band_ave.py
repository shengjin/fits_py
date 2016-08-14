#!/usr/bin/env python2
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

hdulist = fits.open('LkCa15_12CO_cube.image.fits')
hdulist.info()
hdu = hdulist[0]
nq,nc,m,n=hdu.data.shape
hdu.header
print nq, nc, m, n

print hdu.data.min()
print hdu.data.max()

#plt.imshow(hdu.data[0,0,:,:], origin='lower')
#plt.imshow(hdu.data[0,i*2,:,:], origin='lower')
#plt.show()
#hdu.writeto('lat_background_model_slice.fits')
#hdulist.writeto('lat_background_model_slice_allhdus.fits')

hdu.data = hdu.data #+0.0598635
print hdu.data.min()
print hdu.data.max()

sum_channels = np.zeros([m,n], dtype=np.float64)
print sum_channels.shape

RST_freq = 3.457959900000E+11
freq_delt = 2.422247656860E+05
cc = 299792.458   # km/s
dv = freq_delt/RST_freq*cc

pixAU=2.79999863757

for i in range(nc):
    sum_channels = sum_channels + hdu.data[0,i,:,:]*dv

print np.sum(sum_channels, axis=None)

mx = np.arange(1024)
ny = np.arange(1024)
mx = (mx-513)*pixAU
ny = (ny-504)*pixAU

plt.pcolormesh(mx,ny, sum_channels[:,:]*1000, vmin=0, vmax=800)
#plt.ylim( [(-600/pixAU)+504,  (600/pixAU)+504] )
#plt.xlim( [(-600/pixAU)+513,  (600/pixAU)+513] )
plt.xlim( [-600,  600] )
plt.ylim( [-600,  600] )
#plt.colorbar()
plt.ylabel('AU')
plt.xlabel('AU')
plt.colorbar(label='m Janskey/beam')
plt.savefig('mom0.png')




"""
#sizepix in AU
pixAU=2.79999863757

n_row=9
n_col=10

plt, axes = plt.subplots(n_row, n_col, figsize=(15, 12), dpi=80, sharex=True, sharey=True)


for i in range(n_row):
    for j in range(n_col):
        print i,j
        im=axes[i,j].imshow(hdu.data[0,(n_col*i+j),:,:]*1000.0, origin='lower')
        # center y:504 x:513
        axes[i,j].set_ylim(224,784)
        axes[i,j].set_xlim(303,723)
        #im=axes[i,j].pcolormesh(im_x_au, im_y_au, ImageJyppix_scaled[:,:,n_col*i+j], vmin=v_min, vmax=v_max, cmap='jet')
        #im=axes[i,j].pcolormesh(im_x_au, im_y_au, ImageJyppix_scaled[:,:,n_col*i+j], vmin=v_min, vmax=v_max, cmap='RdBu')
        #im=axes[i,j].pcolormesh(im_x_au, im_y_au, (ImageJyppix_scaled[:,:,n_col*i+j]-Image_dust_Jyppix_scaled[:,:,0]), vmin=v_min, vmax=v_max, cmap='RdBu')

        axes[i,j].tick_params(axis='both', which='major', labelsize=8)
        #tick.label.set_fontsize(14) 
        #axes[i,j].axis([im_x_au.min()-axisadd, im_x_au.max()+axisadd, im_y_au.min()-axisadd, im_y_au.max()+axisadd])
        axes[i,j].set_yticks( [(-600/pixAU)+504, (-400/pixAU)+504, (-200/pixAU)+504, 0+504, (200/pixAU)+504, (400/pixAU)+504, (600/pixAU)+504] )
        axes[i,j].set_xticks( [(-600/pixAU)+513, (-400/pixAU)+513, (-200/pixAU)+513, 0+513, (200/pixAU)+513, (400/pixAU)+513, (600/pixAU)+513] )
        axes[i,j].set_xticklabels(['-600', '-400', '-200', '0', '200', '400', '600'],rotation=90)
        axes[i,j].set_yticklabels(['-600', '-400', '-200', '0', '200', '400', '600'])
        if j==0:
            axes[i,j].set_ylabel('AU')
        axes[i,j].set_xlabel('AU')



plt.subplots_adjust(bottom=0.1, right=0.9, top=0.9, hspace=0.0, wspace=0.0)
cax = plt.add_axes([0.92, 0.1, 0.02, 0.8])
plt.colorbar(im, cax=cax, label='m Janskey/beam')

plt.savefig('gas_lines_radmc.png')

plt.clf()

"""










