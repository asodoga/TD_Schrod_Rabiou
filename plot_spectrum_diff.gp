set terminal pngcairo
set output 'spectrum_stdnb70_Hnb5.png'
set title "spectrum_diff"
set xlabel "temps"
set ylabel "spectrum"
set boxwidth 2
 plot 'file_spectrum_stdnb70.txt'  u 1:3 w l lw 2 t 'spectrum_stdnb70 op=3',\
 'file_spectrum_Hnb5.txt'  u 1:3 w p pt 7 ps 0.6 t 'spectrum_Hnb5 op=3'
   