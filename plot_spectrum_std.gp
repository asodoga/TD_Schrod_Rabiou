set terminal pngcairo
set output 'spectrum_std.png'
set title "spectrum_std"
set xlabel "temps"
set ylabel "spectrum"
set boxwidth 2
 plot 'file_spectrum_stdnb50.txt'  u 1:3 w l lw 2 t 'spectrum_stdnb50 op=3'
   