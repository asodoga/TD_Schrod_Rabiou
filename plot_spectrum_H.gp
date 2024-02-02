set terminal pngcairo
set output 'spectrum_H.png'
set title "spectrum_H"
set xlabel "temps"
set ylabel "spectrum"
set boxwidth 2
 plot 'file_spectrum_Hnb50.txt'  u 1:3 w l lw 2 t 'spectrum_Hnb50 op=3'
   