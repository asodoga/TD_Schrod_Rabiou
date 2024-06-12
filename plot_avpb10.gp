set term pdf size 40,30 font "Arial, 90"
set output 'avpnb10.pdf'
 

set ylabel "{P_t}(ua)" font "Arial,90"
set xlabel "Time(ua)" font "Arial,90"

set border linewidth 18.0
set boxwidth  3
set style fill solid 0.5 noborder
plot [-0.5:60][-2.5:3] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/imp_k_hagedorn_taylor.txt' every 1 u 1:2 w lp lw 10 pt 7   ps 3 t "Pt_{1}" ,\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/imp_k_hagedorn_taylor.txt'  every 10 u 1:3 w lp lw 10  pt 5  ps 3 t "Pt_{2}"

   
