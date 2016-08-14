#!/usr/bin/env python2
from astropy.io import fits
import matplotlib.pyplot as plt

import numpy as np

import matplotlib as mpl
mpl.use('Agg')  # to avoid window

from matplotlib.colors import LogNorm

from matplotlib import cm
import math
from math import sin, cos

# LkCa 15 continuum
noise_beam_flux = 9.6e-5 # Jy/beam

# C12O
#noise_beam_flux = 2.1e-2 # Jy/beam km/s
# C13O
#noise_beam_flux = 1.1e-2 # Jy/beam km/s
# 18CO
#noise_beam_flux = 1.1e-2 # Jy/beam km/s

# LkCa 15 CO

# radius_wanted
n_rad = 60
dr = 10
inneravoid = 0
gap_width = 10.0 # AU

n_image = 0

# inclination of the disk
incli = 51.0 #incli = 46.72

# position angle
# use it to scale the center Radius VS T_B
posang = 151.0 #posang = 138.02

# rotation angle used in the center ellipse to
rotang = posang
#rotang = 90-posang

# center of the ellipse in pixel
e_h = 513.0
e_k = 504.0

dpc = 140.0  # distance in pc
n_image = 0  # the n_th image in image.out

#debug = True
debug = False

AU = 1.496e13
pc = 3.08572e18 # pc in cm


#hdulist = fits.open('LkCa15_13CO_cube.image.mom0.fits')
#hdulist = fits.open('LkCa15_C18O_cube.image.mom0.fits')
hdulist = fits.open('LkCa15.continuum_map_r0.5.image.fits')
print "hdulist.info"
hdulist.info()
print ""

hdu = hdulist[0]
m,n,o,p=hdu.data.shape
print "hdu.data.shape"
print m,n,o,p
print ""

gaus_ma = hdu.header['BMAJ']
gaus_mi = hdu.header['BMIN']

print 'NAXIS'
print hdu.header['NAXIS']
print 'NAXIS1'
print hdu.header['NAXIS1']
print 'NAXIS2'
print hdu.header['NAXIS2']
print 'NAXIS3'
print hdu.header['NAXIS3']
print 'NAXIS4'
print hdu.header['NAXIS4']
print ""
print 'BMAJ'
print hdu.header['BMAJ']
print 'BMIN'
print hdu.header['BMIN']
print 'BPA'
print hdu.header['BPA']
print ""
print ['BUNIT']
print hdu.header['BUNIT']
print ""
print 'CTYPE1'
print hdu.header['CTYPE1']
print 'CDELT1'
print hdu.header['CDELT1']
print ""
print 'CTYPE2'
print hdu.header['CTYPE2']
print 'CDELT2'
print hdu.header['CDELT2']
print ""
print 'CTYPE3'
print hdu.header['CTYPE3']
print 'CRVAL3'
print hdu.header['CRVAL3']
print ""
print 'CTYPE4'
print hdu.header['CTYPE4']
print ""



# Calculate the effective beam area in px given the Gaussian FWHMs      
degpix = hdu.header['CDELT1']
degpiy = hdu.header['CDELT2']
# Convert FWHMs to sigma
gfactor = 2.0*math.sqrt(2*math.log(2))  # FWHM=gfactor*sigma
sigma_maj = gaus_ma/gfactor
sigma_min = gaus_mi/gfactor
# Calculate the beam area in the input units and in pixels
beam_area_sq = 2.0*math.pi*sigma_maj*sigma_min
pix_area_sq = abs(degpix*degpiy)
beam_area_pix = beam_area_sq/pix_area_sq
print "beam_area_pix    "
print beam_area_pix    


nx = hdu.header['NAXIS1']

pixel_in_rad = abs(degpix)*math.pi/180.0
pixely_in_rad = abs(degpiy)*math.pi/180.0
sizepix = math.sin(pixel_in_rad)*dpc*pc
sizepiy = math.sin(pixely_in_rad)*dpc*pc

print "nx, sizepix"
print nx, sizepix
print "sizepix in AU"
print sizepix/AU
print sizepiy/AU



##########################################
# define azimuthal extract function
#  could be ellipse or circle
##########################################

def azimuthal_Jy_avg(image, gap_au, gap_width, e_h, e_k, incli, rotang, sizepix):
    # input: image, gap_au, gap_width, e_h, e_k, incli, rotang, sizepix, nx
    # output: fin_vis
    ##################################
    # screen out the point in the ring at [i_in, i_out]

    au = 1.496e13

    fin0_x = []
    fin0_y = []
    fin0_vis = []
    
    # parameters for the center ellipse
    gap_min = gap_au - gap_width*0.5
    gap_max = gap_au + gap_width*0.5
    if gap_min < 0.0:
        gap_min = 0.0
    
    # assuming sizepix_x = sizepix_y
    e_a_grid_max = gap_max*au/sizepix  # long semi-axis
    e_a_grid_min = gap_min*au/sizepix  # long semi-axis
    
    inclination = math.radians(incli)
    e_b_grid_max = e_a_grid_max * cos(inclination)     # short semi-axis
    e_b_grid_min = e_a_grid_min * cos(inclination)     # short semi-axis

    rotang = math.radians(rotang)

    m = image.shape[0]
    # convert integer to float in order to make sure we find
    #    the center of the image.
    # image.out: 1) 20 grids
    #               python array: 0, 1, ..., 19
    #               center is at point 10
    #            2) 19 grids
    #               python array: 0, 1, ..., 18
    #               center is at point 9.5
    i_2_au = sizepix/au

    # NOTE: i,j change position
    for ii in range(m):
        j=float(ii)
        for jj in range(m):
            i=float(jj)
            if (e_a_grid_min == 0.0) and (e_b_grid_min == 0.0):
                if ( ((i-e_h)*cos(rotang)+(j-e_k)*sin(rotang))**2/e_a_grid_max**2 + ((i-e_h)*sin(rotang)-( j-e_k)*cos(rotang))**2/e_b_grid_max**2 <= 1.0**2):
                    fin0_x.append(i)
                    fin0_y.append(j)
                    fin0_vis.append(image[ii,jj])
            else:
                if ( ((i-e_h)*cos(rotang)+(j-e_k)*sin(rotang))**2/e_a_grid_max**2 + ((i-e_h)*sin(rotang)-( j-e_k)*cos(rotang))**2/e_b_grid_max**2 <= 1.0**2) and ( ((i-e_h)*cos(rotang)+(j-e_k)*sin(rotang))**2/e_a_grid_min**2 + ((i-e_h)*sin(rotang)-(j-e_k)*cos(rotang))**2/e_b_grid_min**2 > 1.0**2) :
                    fin0_x.append(i)
                    fin0_y.append(j)
                    fin0_vis.append(image[ii,jj])
                #print nhalfpix, i, j
    fin_x = np.asarray(fin0_x)
    fin_y = np.asarray(fin0_y)
    fin_vis = np.asarray(fin0_vis)
    n_fin = fin_x.shape[0]
    if n_fin > 0: 
        np.savetxt('x_y_vis', np.transpose([fin_x,fin_y,fin_vis]))
    #total=np.sum(fin_vis), "\n"
    #avg = total/n_fin
    return fin_vis

##########################################
##########################################


RST_freq = 3.457959900000E+11
freq_delt = 2.422247656860E+05
cc = 299792.458   # km/s
dv = freq_delt/RST_freq*cc

Image2D = np.zeros([nx,nx], dtype=np.float64)
Image2D = hdu.data[0,0,:,:]
#Image2D = Image2D+0.0598635*90*dv


print Image2D[e_h, e_k]

r_Jy = np.zeros([n_rad,3], dtype=np.float64)

for i in range(n_rad):
    gap_au = float(i)*dr+inneravoid
    r_Jy[i,0] = gap_au
    avg = azimuthal_Jy_avg(Image2D, gap_au, gap_width, e_h, e_k, incli, rotang, sizepix)
    
    # parameters for the center ellipse
    gap_min = gap_au - gap_width*0.5
    gap_max = gap_au + gap_width*0.5
    if gap_min < 0.0:
        gap_min = 0.0
    
    # assuming sizepix_x = sizepix_y
    e_a_grid_max = gap_max*AU/sizepix  # long semi-axis
    e_a_grid_min = gap_min*AU/sizepix  # long semi-axis
    inclination = math.radians(incli)
    e_b_grid_max = e_a_grid_max * cos(inclination)     # short semi-axis
    e_b_grid_min = e_a_grid_min * cos(inclination)     # short semi-axis
    area_max_pix = math.pi*e_a_grid_max*e_b_grid_max 
    area_min_pix = math.pi*e_a_grid_min*e_b_grid_min 
    area_ring_pix = area_max_pix - area_min_pix
    print area_max_pix , area_min_pix
    print area_ring_pix,beam_area_pix
    print noise_beam_flux, math.sqrt(area_ring_pix/beam_area_pix)
    r_Jy[i,2] = noise_beam_flux/(math.sqrt(area_ring_pix/beam_area_pix))
    
    n_num = avg.shape[0]
    if n_num == 0:
        r_Jy[i,1] = 0
        print "WARNNING: no points found around ", r_Jy[i,0], " AU!!"
    else:
        f_n_num = float(n_num)
        #print avg
        total = np.sum(avg)
        #print total
        r_Jy[i,1] = total/f_n_num
        print r_Jy[i,0], r_Jy[i,1], r_Jy[i,2]
    
    

#######################################################
#####  Make some Plot 
#########################################


np.savetxt('AU_fluxJy_noissJy', np.transpose([r_Jy[:,0], r_Jy[:,1], r_Jy[:,2]]))

plt.xlabel('AU')
plt.ylabel("m Janskey / beam")
#plt.ylim(-180,180)
#plt.xlim(0,n_rad+inneravoid)
plt.errorbar(r_Jy[:,0], r_Jy[:,1]*1000, xerr=0, yerr=r_Jy[:,2]*1000)
#plt.pcolormesh(im_x_au, im_y_au, ImageJyppix_scaled*1000, vmin=1e-6, vmax=0.002,norm=LogNorm(),cmap='RdBu')
plt.savefig('azimuthalavg_fluxJy.png')
plt.clf()

pixAU=2.79999863757

mx = np.arange(1024)
ny = np.arange(1024)
mx = (mx-513)*pixAU
ny = (ny-504)*pixAU

plt.pcolormesh(mx,ny,hdu.data[0,0,:,:]*1000, vmin=0, vmax=18) #, origin='lower')
plt.xlim(-100,100)
plt.ylim(-100,100)
#plt.imshow(hdu.data[0,i*2,:,:], origin='lower')
#plt.show()
#hdu.writeto('lat_background_model_slice.fits')

plt.ylabel('AU')
plt.xlabel('AU')
plt.colorbar(label='m Janskey/beam')
plt.savefig('dust.png')















