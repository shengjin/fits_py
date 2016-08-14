#!/usr/bin/env python2

import math

AU = 1.496e13
#pixel_in_deg = 5.555555555556E-06    # CUNIT2  = 'deg
pixel_in_deg = 9.722222222222E-07
pixel_in_rad = pixel_in_deg*math.pi/180.0
dpc = 140.0     # distance of pc
pc = 3.08572e18 # pc in cm
sizeAU = math.sin(pixel_in_rad)*dpc*pc/AU

print pixel_in_deg
print pixel_in_rad
print pixel_in_deg*3600
print pixel_in_rad*206264.806247
print math.sin(pixel_in_rad)
print math.sin(pixel_in_rad)*dpc*pc
print sizeAU
