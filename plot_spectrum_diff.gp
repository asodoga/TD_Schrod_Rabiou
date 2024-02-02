set terminal pngcairo
set output 'spectrum_diff.png'
set title "spectrum_diff"
set xlabel "temps"
set ylabel "spectrum"
set boxwidth 2
 plot 'file_spectrum_stdnb50.txt'  u 1:3 w l lw 2 t 'spectrum_stdnb50 op=3',\
 'file_spectrum_Hnb50.txt'  u 1:3 w p pt 7 ps 0.6 t 'spectrum_Hnb50 op=3'
   