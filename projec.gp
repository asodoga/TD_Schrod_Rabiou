set terminal pngcairo
set output 'projection.png'
set title "projection"
set xlabel "Time"
set ylabel "projection"
set boxwidth 2
#set logscale y 2
#set key left bottom

 plot [:][:]'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb100/proj'  w l lw 2 t 'projection  nb=10' 