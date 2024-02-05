set terminal pngcairo
set output 'spectrum_diff_stdnb70_Hnb45.png'
set title "spectrum_diff"
set xlabel "temps"
set ylabel "spectrum"
set boxwidth 2
 plot 'file_spectrum_stdnb70_Hnb45.txt'  w l lw 4 t 'spectrum_diff_stdnb70_Hnb45 op=3'