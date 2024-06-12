reset
set term pdf size 40,30 font "Arial, 90"

set output 'avqnb10.pdf'
set border linewidth 18.0

set ylabel "{Q_t} (ua)" font "Arial,90"
set xlabel "Time (ua)" font "Arial,90"


plot [-0.5:60][-2.5:3] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/Qt_hagedorn_taylor.txt' every 1 u 1:2 w lp lw 10 pt 7 ps 3 t"{Qt_{1}}" ,\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb10/Qt_hagedorn_taylor.txt' every 10 u 1:3 w lp lw 10 pt 5 ps 3  t"{Qt_{2}}"


   
