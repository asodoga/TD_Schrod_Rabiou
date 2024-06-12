set term pdf size 35,25 font "Arial, 90"
set output 'Harmonic-N13-error.pdf'

set border linewidth 18.0

set ylabel "{/Symbol D}{N_{1,3}}(t)" font "Arial,90"
set xlabel "Time (ua)" font "Arial,90"
set logscale y
set boxwidth  3

plot [5:60][:0.2] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5/Norm_13_hagedorn_taylor.txt' u 1:2 w l lw 10 t " nb=5",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/Norm_13_hagedorn_taylor.txt' every 25 u 1:2 w lp pt 7  lw 10 ps 7 t " nb=10 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb20/Norm_13_hagedorn_taylor.txt' every 15 u 1:2 w lp pt 5  lw 10 ps 7 t " nb=20 ",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb30/Norm_13_hagedorn_taylor.txt' every 10 u 1:2 w lp pt 3  lw 10 ps 7 t " nb=30 "

#pause -1
   
