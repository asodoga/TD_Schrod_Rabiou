set term pdf size 40,30 font "Arial, 90"
set output 'Harmonic-r-pics-error.pdf'
set tics font "Arial,90"
 
set border linewidth 18.0

set ylabel "{R_C}(t)" font "Arial,90"
set xlabel "Time (ua)"font "Arial,90"
set logscale y
set boxwidth  3

plot[5:60][:4] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5/file_norm_pics_hagedorn_taylor.txt' u 1:3 w l lw 10 t"nb=5",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/file_norm_pics_hagedorn_taylor.txt' every 10 u 1:3 w lp pt 7  lw 10 ps 7 t"nb=10",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb20/file_norm_pics_hagedorn_taylor.txt' every 10 u 1:3 w lp pt 5  lw 10 ps 7 t"nb=20",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb30/file_norm_pics_hagedorn_taylor.txt' every 5 u 1:3 w lp pt 3   lw 10 ps 7 t"nb=30"


#pause -1
   
