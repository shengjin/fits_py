import numpy as np
import aplpy

# Convert all images to common projection
aplpy.make_rgb_cube(['m1.fits', 'i3.fits', 'i2.fits'], 'rgb.fits')

# Make 3-color image
aplpy.make_rgb_image('rgb.fits', 'rgb.png',
                     vmin_r=20, vmax_r=400,
                     vmin_g=0, vmax_g=150,
                     vmin_b=-2,vmax_b=50,
                     embed_avm_tags=True)

# Make a plot similarly to before
fig = aplpy.FITSFigure('rgb.png')
fig.show_rgb()
fig.show_contour('sc.fits', cmap='jet', levels=np.linspace(0.0, 1.0, 10))
fig.ticks.set_color('white')
fig.tick_labels.set_font(size='small')
fig.add_grid()
fig.grid.set_alpha(0.5)
fig.save('aplpy_rgb_plot.png')

