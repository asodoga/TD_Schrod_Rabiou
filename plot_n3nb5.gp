set term pdf size 40,30 font "Arial, 90"
set output 'Harmonic-N3nb5-error.pdf'
 
set border linewidth 18.0

set ylabel "{/Symbol D}{N_{3}}(t)" font "Arial,90"
set xlabel "Time (ua)" font "Arial,90"
set logscale y
set boxwidth  3

plot[8:60][:0.2] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.4/Norm_hagedorn_taylor.txt' u 1:(1-$2) w l lw 10 t"dt=0.4" ,\
               '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.3/Norm_hagedorn_taylor.txt' every 10 u 1:(1-$2) w lp pt 7 lw 10 ps 7 t"dt=0.3",\
               '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.2/Norm_hagedorn_taylor.txt' every 20 u 1:(1-$2) w lp pt 5 lw 10 ps 7 t"dt=0.2",\
               '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.1/Norm_hagedorn_taylor.txt' every 50 u 1:(1-$2) w lp pt 3 lw 10 ps 7 t"dt=0.1",\
               '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb5dt0.01/Norm_hagedorn_taylor.txt' every 100 u 1:(1-$2) w lp pt 12 lw 10 ps 7 t"dt=0.01"

#pause -1
   
