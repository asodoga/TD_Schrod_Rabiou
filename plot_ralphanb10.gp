set term pdf size 40,30 font "Arial, 90"
set output 'ralpha-plotnb10.pdf'
set tics font "Arial,90"
 
set border linewidth 18.0

set ylabel "{Re({/Symbol a_t})}" font "Arial,90"
set xlabel "Time (ua)" font "Arial,90"
set boxwidth  3

plot [-1 :60][0.5:1.6] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/alpha_hagedorn_taylor.txt' every 1 u 1:2 w lp lw 7 pt 7  ps 3 t "{Re({/Symbol a}_1)}" ,\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/alpha_hagedorn_taylor.txt' every 10 u 1:3 w lp lw 7 pt 5  ps 3  t "{Re({/Symbol a}_2)}"


   
