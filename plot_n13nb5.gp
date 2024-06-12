set term pdf size 35,25 font "Arial, 90"
set output 'Harmonic-N13nb5-error.pdf'
set tics font "Arial,90"
 
set border linewidth 18.0

set ylabel "{/Symbol D}{N_{1,3}}(t)" font "Arial,90"
set xlabel "Time (ua)" font "Arial,90"
set logscale y
set boxwidth  3

plot [5:60][:0.01] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.4/Norm_13_hagedorn_taylor.txt' u 1:2 w l lw 10 t"dt=0.4",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.2/Norm_13_hagedorn_taylor.txt' every 20 u 1:2 w lp pt 7 lw 10  ps 7 t"dt=0.2",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.1/Norm_13_hagedorn_taylor.txt' every 20 u 1:2 w lp pt 5 lw 10 ps 7 t"dt=0.1",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.01/Norm_13_hagedorn_taylor.txt' every 500 u 1:2 w lp pt 3 lw 10 ps 7 t"dt=0.01"


#pause -1
   
