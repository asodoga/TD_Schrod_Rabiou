 
 set terminal pngcairo
set output 'Imp_k_diff.png'
 set title "Impk diff"
set xlabel "temps"
set ylabel "Impk"
set boxwidth 2
plot "Impk_diff.txt" u 1:($2-$5) w l lw 4  t 'Imp_k diff '
