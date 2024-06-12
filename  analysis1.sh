#!/bin/bash
#./calc_spectrum_std.sh 70

# cat >  plot_spectrum_stdnb70stdnb5Hnb5.gp<< EOF
cat >  plot_spectrum_stdnb70stdnb5Hnb5.gp<< EOF

set terminal pngcairo
set output 'spectrum_stdnb70stdnb5Hnb5.png'
set title "spectrum_std"
set xlabel "E[ua]"
set ylabel "spectrum"
set boxwidth 2
 plot 'file_spectrum_stdnb70.txt'  u 1:3 w l lw 2 t 'spectrum_stdnb70 op=3',\
 plot 'file_spectrum_stdnb5.txt'  u 1:(-$3) w l lw 2 t 'spectrum_stdnb5 op=3',\
 plot 'file_spectrum_Hnb5.txt'  u 1:(-$3) w l lw 2 t 'spectrum_Hnb5 op=3'
   
EOF








