# coding: utf-8
import aplpy
f = aplpy.FITSFigure('LkCa15.continuum_map_r0.5.image.fits')
f.show_grayscale()
f.show_contour('LkCa15.continuum_map_r0.5.image.fits', levels=10)
f.add_grid()
f.add_scalebar(0.0003, '0.5pc', color='white')
f.save('my_first_plot.eps')
