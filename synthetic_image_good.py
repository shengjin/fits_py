import pyfits
import numpy as np
from scipy.ndimage import gaussian_filter

# Create empty image
nx, ny = 512, 512
image = np.zeros((ny, nx))

# Set number of stars
n = 10000

# Generate random positions
r = np.random.random(n) * nx
theta = np.random.uniform(0., 2. * np.pi, n)

# Generate random fluxes
f = np.random.random(n) ** 2

# Compute position
x = nx / 2 + r * np.cos(theta)
y = ny / 2 + r * np.sin(theta)

# Add stars to image
# ==> First for loop and if statement <==
for i in range(n):
    if x[i] >= 0 and x[i] < nx and y[i] >= 0 and y[i] < ny:
        image[y[i], x[i]] += f[i]

# Convolve with a gaussian
image = gaussian_filter(image, 1)

# Add noise
image += np.random.normal(3., 0.01, image.shape)

# Write out to FITS image
pyfits.writeto('cluster.fits', image, clobber=True)
