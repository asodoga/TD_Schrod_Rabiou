set term pdf size 40,30 font "Arial, 90"
set output 'Harmonic-r-pics-errornb5.pdf'
set tics font "Arial,90"
 set border linewidth 18.0

set ylabel "{R_C}(t) " font "Arial,90"
set xlabel "Time (ua)" font "Arial,90"
set logscale y
set boxwidth  3


plot[5:60][:4] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.4/file_norm_pics_hagedorn_taylor.txt' u 1:3 w l lw 10 t"dt=0.4",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.3/file_norm_pics_hagedorn_taylor.txt' every 10 u 1:3 w lp pt 7 lw 10 ps 7 t"dt=0.3",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.2/file_norm_pics_hagedorn_taylor.txt' every 10 u 1:3 w lp pt 5 lw 10 ps 7 t"dt=0.2",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.1/file_norm_pics_hagedorn_taylor.txt' every 50 u 1:3 w lp pt 3 lw 10 ps 7 t"dt=0.1",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.01/file_norm_pics_hagedorn_taylor.txt' every 100 u 1:3 w lp pt 12 lw 10  ps 7 t"dt=0.01"

   
