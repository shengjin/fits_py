set terminal postscript eps enhanced defaultplex \
   leveldefault color colortext \
   dashed dashlength 1.0 linewidth 2.0 butt \
   palfuncparam 2000,0.003 \
   "SFRM1000" 28 fontfile "/home/jin/fonts/sfrm1000.pfb" 
set output "test.eps"


set bmargin 0
set rmargin 0
set lmargin 0
set tmargin 0



# linetype color for label
# heavy color
set style line 1  linetype 1 linecolor rgb "#000000"  linewidth 1.500 pointtype  3 pointsize default #blue
#set style line 1  linetype 1 linecolor rgb "#1d4599"  linewidth 1.500 pointtype  3 pointsize default #blue
set style line 2  linetype 3 linecolor rgb "#1d4599"  linewidth 1.2 pointtype  3 pointsize default #blue
#set style line 2  linetype 3 linecolor rgb "#1d4599"  linewidth 0.5 pointtype  3 pointsize default #blue
set style line 5  linetype 1 linecolor rgb "#000000"  linewidth 1.000 pointtype  3 pointsize default #green
#set style line 5  linetype 1 linecolor rgb "#11ad34"  linewidth 1.000 pointtype  3 pointsize default #green
set style line 3  linetype 5 linecolor rgb "#000000"  linewidth 1.000 pointtype  3 pointsize default dashtype '--'
set style line 4  linetype 1 linecolor rgb "#e69f17"  linewidth 1.500 pointtype  3 pointsize default #orange

xtotalsize = 1.0
ytotalsize = 1.5
xsizeadj=0.00

set size xtotalsize-xsizeadj,ytotalsize

#set size ratio 0.65


set multiplot 


ysize =  0.8
xsize =  0.8
set size xsize,ysize

#set xrange [ 0.00000 : 13.0000 ] noreverse nowriteback

set origin 0.16,0.61
set format x ""
set format y "%.1f"
#set format y "10^{%L}"
set mytics 2
set yrange [0:12]
set xrange [-0.6:0.6]
set label 1 "run1" at graph 0.18, graph 0.90 right font "Helvetica,28" textcolor lt -1
set ylabel "Flux density [ mJy/beam ]" offset 1.0,-3.5 font ",28"
plot '../image_observation/radial_profile.ascii' u (($1-292)*0.0021428877):(($2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12)/11*1000) w l ls 3 t "Observation" \
, 'diffmass_RmIm' u (($1-292)*0.0021428877):(($2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12)/11*1000) w l ls 1 t "Model"


unset title
ysize =  0.3
xsize =  0.8
set size xsize,ysize

set origin 0.16,0.25
set format x "%.1f"
set format y "%.1f"
set yrange [-1:1]
set xrange [-0.6:0.6]
#set ytics 0,1,6
set mytics 1
unset label
set xlabel "Radius [ arcsec ]" font ",28" offset 0,-0.5
unset ylabel
set key font ",20"
#set key at 12.2,0.10 font ",17"
plot 'diffmass_RdId' u (($1-292)*0.0021428877):(($2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12)/11*1000) w l ls 5 t "Obs-Mod" 


unset multiplot

