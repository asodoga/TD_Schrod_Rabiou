set terminal pngcairo
#set terminal postscript enhanced
set output 'Harmonic-pics-error.png'
set tics font "Times-news-romain,16"
 
set origin 0,0
#set size 1.5,1.75
set multiplot
set border linewidth 2.0

set ylabel "Intensity (ua)" font "Times-news-romain,14"
set xlabel "Time (ua)" font "Times-news-romain,14"
set title "  error parameter with {nb_1}={nb_2} = 5 " font "Times-news-romain,14"
set xtics 6

set tics out nomirror
set logscale y
set boxwidth  3

plot [-3:63][:] '/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb25/file_norm_pics_hagedorn_taylor.txt' u 1:2 w l lw 2 t "{/Symbol D}{C_1}(t)",\
'/Users/issa/Developpements/TD_Schrod_Rabiou/results_H_nb25/file_norm_pics_hagedorn_taylor.txt' u 1:3 w l lw 2 t "{R_C}(t)"







#pause -1