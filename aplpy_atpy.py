import numpy as np
import aplpy
import atpy

# Create a new figure
fig = aplpy.FITSFigure('i3.fits')

# Show the colorscale
fig.show_colorscale()

# Add contours
fig.show_contour('sc.fits', cmap='jet', levels=np.linspace(0.0, 1.0, 10))

# Make ticks white
fig.ticks.set_color('white')

# Make labels smaller
fig.tick_labels.set_font(size='small')

# Overlay a grid
fig.add_grid()
fig.grid.set_alpha(0.5)

# Add a colorbar
fig.add_colorbar()

# Use ATpy to read in an IRSA table
tab = atpy.Table('2mass.tbl')
tab_bright = tab.where(tab['j_m'] < 13.)

# Plot markers
fig.show_markers(tab['ra'], tab['dec'], marker='+', edgecolor='white')
fig.show_markers(tab_bright['ra'], tab_bright['dec'], marker='o', edgecolor='white')

# Save image for publication
fig.save('aplpy_plot.png')
