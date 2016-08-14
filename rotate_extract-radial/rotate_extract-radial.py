ia.open('HLTau_B7cont_av20_mscale_ap_shifted.image.fits')
imr=ia.rotate(pa='-41.98deg')
ok = imr.tofits('HLTau_B7cont_av20_mscale_ap_shifted.image_rot42.fits', overwrite=true)
imr.done()
ia.close()

val = imval('HLTau_B7cont_av20_mscale_ap_shifted.image_rot42.fits', box='730,1019,1318,1029', stokes='I')
f = open('radial_profile.ascii', 'w')

for m in range (0, 589):
    
    av = 0.0
    s = str(m)+' '
    f.write(s)

    for n in range (0, 10):
        av = av+val['data'][m][n]
        s = str(val['data'][m][n])+' '
        f.write(s)
    
    av = av/10.
    s = str(av)
    f.write(s)
    f.write('\n')
    
f.close()
